library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

-- Show-ahead (first-word fall-through) FIFO: o_rd_data always presents the
-- oldest word (flagged by o_rd_valid); i_rd_en pops it. A word written into
-- an empty FIFO reaches o_rd_data one cycle later (plain register semantics).
--
-- Structure (oldest -> youngest): head -> mid -> rd_q -> mem.
--   * mem + rd_q form a plain synchronous-read RAM: rd_q is loaded from
--     fifo_mem and from NOTHING else, so Vivado maps it onto a block-RAM
--     primitive's output register (a mux in front of that register would
--     force distributed RAM).
--   * head/mid are two fabric registers providing the show-ahead interface
--     and absorbing the 2-cycle mem->rd_q->head refill latency, so streaming
--     pops run at 1 word/cycle with no bubbles.
--   * writes append at the youngest occupied end; they bypass directly into
--     head/mid only when rd_q and mem are empty (order preserved), which
--     keeps the 1-cycle write-to-head behavior a depth-1 delay line needs.
entity fifo_fwft is
	generic (
		G_DATA_WIDTH : integer := 36;
		G_DEPTH      : integer := 1024;
		G_RAM_STYLE  : string := "auto"
	);
	port (
		i_clk 		: in std_logic;
		i_reset 	: in std_logic;

		i_wr_en 	: in std_logic;
		i_wr_data : in std_logic_vector(G_DATA_WIDTH-1 downto 0);

		i_rd_en 	 : in std_logic;   -- pop: consume the presented head word
		o_rd_data  : out std_logic_vector(G_DATA_WIDTH-1 downto 0);
		o_rd_valid : out std_logic;

		o_full 		: out std_logic
	);
end entity fifo_fwft;

architecture rtl of fifo_fwft is

	-- ceil(log2(n)), minimum of 1 bit.
	function clogb2(n : integer) return integer is
		variable res : integer := 0;
		variable v   : integer := n - 1;
	begin
		while v > 0 loop
			res := res + 1;
			v   := v / 2;
		end loop;
		if res = 0 then
			res := 1;
		end if;
		return res;
	end function;

	constant C_ADDR_W : integer := clogb2(G_DEPTH);      -- address bits
	constant C_CNT_W  : integer := clogb2(G_DEPTH + 1);  -- occupancy 0..G_DEPTH


	type t_mem is array (0 to G_DEPTH-1) of std_logic_vector(G_DATA_WIDTH-1 downto 0);
	signal fifo_mem : t_mem;

	attribute ram_style : string;
	attribute ram_style of fifo_mem : signal is G_RAM_STYLE;

	-- Pointers wrap at G_DEPTH (not necessarily a power of two); mem_cnt is
	-- the number of words in the memory (head/mid/rd_q are extra capacity).
	signal wr_ptr  : unsigned(C_ADDR_W-1 downto 0) := (others => '0');
	signal rd_ptr  : unsigned(C_ADDR_W-1 downto 0) := (others => '0');
	signal mem_cnt : unsigned(C_CNT_W-1 downto 0)  := (others => '0');

	signal rd_q : std_logic_vector(G_DATA_WIDTH-1 downto 0);  -- BRAM output reg
	signal rq_v : std_logic := '0';

	signal head, mid     : std_logic_vector(G_DATA_WIDTH-1 downto 0);
	signal head_v, mid_v : std_logic := '0';

	signal full_s : std_logic;

	-- combinational schedule for this cycle
	signal pop         : std_logic;
	signal head_free   : std_logic;   -- head slot open after this edge
	signal mid_to_head : std_logic;
	signal rq_to_head  : std_logic;
	signal rq_to_mid   : std_logic;
	signal mem_rd      : std_logic;
	signal tail_empty  : std_logic;   -- no older word in rd_q or mem
	signal wr_to_head  : std_logic;
	signal wr_to_mid   : std_logic;
	signal wr_to_mem   : std_logic;

begin

	full_s <= '1' when mem_cnt = G_DEPTH else '0';
	pop    <= i_rd_en and head_v;

	head_free   <= (not head_v) or pop;
	mid_to_head <= mid_v and head_free;
	rq_to_head  <= rq_v and head_free and not mid_v;
	rq_to_mid   <= rq_v and not rq_to_head and ((not mid_v) or mid_to_head);

	-- prefetch: keep rd_q loaded whenever the memory holds data and rd_q is
	-- (or is becoming) free
	mem_rd <= '1' when mem_cnt /= 0 and
	                   (rq_v = '0' or rq_to_head = '1' or rq_to_mid = '1')
	          else '0';

	tail_empty <= '1' when mem_cnt = 0 and rq_v = '0' else '0';

	-- write routing: bypass into the fabric registers only when no older word
	-- sits behind them (order), else append into the memory
	wr_to_head <= i_wr_en and tail_empty and (not mid_v) and head_free;
	wr_to_mid  <= i_wr_en and tail_empty and not wr_to_head and
	              ( ((not mid_v) and head_v and not pop) or mid_to_head );
	wr_to_mem  <= i_wr_en and not wr_to_head and not wr_to_mid and not full_s;

	-- ---- BRAM core: write port + sync-read into the dedicated rd_q ----
	PROC_MEM: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if wr_to_mem = '1' then
				fifo_mem(to_integer(wr_ptr)) <= i_wr_data;
			end if;
			if mem_rd = '1' then
				rd_q <= fifo_mem(to_integer(rd_ptr));
			end if;
		end if;
	end process PROC_MEM;

	-- ---- fabric control: pointers, occupancy, show-ahead registers ----
	PROC_CTRL: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				wr_ptr  <= (others => '0');
				rd_ptr  <= (others => '0');
				mem_cnt <= (others => '0');
				rq_v    <= '0';
				head_v  <= '0';
				mid_v   <= '0';
			else
				if mid_to_head = '1' then
					head <= mid;
					head_v <= '1';
				elsif rq_to_head = '1' then
					head <= rd_q;
					head_v <= '1';
				elsif wr_to_head = '1' then
					head <= i_wr_data;
					head_v <= '1';
				elsif pop = '1' then
					head_v <= '0';
				end if;

				if rq_to_mid = '1' then
					mid <= rd_q;
					mid_v <= '1';
				elsif wr_to_mid = '1' then
					mid <= i_wr_data;
					mid_v <= '1';
				elsif mid_to_head = '1' then
					mid_v <= '0';
				end if;

				rq_v <= mem_rd or (rq_v and not (rq_to_head or rq_to_mid));

				if wr_to_mem = '1' then
					if wr_ptr = G_DEPTH-1 then
						wr_ptr <= (others => '0');
					else
						wr_ptr <= wr_ptr + 1;
					end if;
				end if;
				if mem_rd = '1' then
					if rd_ptr = G_DEPTH-1 then
						rd_ptr <= (others => '0');
					else
						rd_ptr <= rd_ptr + 1;
					end if;
				end if;
				if wr_to_mem = '1' and mem_rd = '0' then
					mem_cnt <= mem_cnt + 1;
				elsif mem_rd = '1' and wr_to_mem = '0' then
					mem_cnt <= mem_cnt - 1;
				end if;
			end if;
		end if;
	end process PROC_CTRL;

	o_rd_data  <= head;
	o_rd_valid <= head_v;
	o_full     <= full_s;

end architecture rtl;
