library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity fifo is
	generic (
		G_DATA_WIDTH : integer := 36;
		G_DEPTH 		 : integer := 1024
	);
	port (
		i_clk 		: in std_logic;
		i_reset 	: in std_logic;

		-- runtime depth: full trips at this occupancy
		i_virtual_depth : in integer range 1 to G_DEPTH := G_DEPTH;

		i_wr_en 	: in std_logic;
		i_wr_data : in std_logic_vector(G_DATA_WIDTH-1 downto 0);

		i_rd_en 	: in std_logic;
		o_rd_data : out std_logic_vector(G_DATA_WIDTH-1 downto 0);

		o_full 		: out std_logic;
		o_almost_full : out std_logic;  -- holds i_virtual_depth-1 samples
		o_empty 	: out std_logic
	);
end entity fifo;

architecture rtl of fifo is

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

	-- Type definitions
	type t_mem is array (0 to G_DEPTH-1) of std_logic_vector(G_DATA_WIDTH-1 downto 0);
	signal fifo_mem : t_mem;

	-- pointers wrap at G_DEPTH, count tells full from empty
	signal wr_ptr : unsigned(C_ADDR_W-1 downto 0) := (others => '0');
	signal rd_ptr : unsigned(C_ADDR_W-1 downto 0) := (others => '0');
	signal count  : unsigned(C_CNT_W-1 downto 0)  := (others => '0');

	signal full_s  : std_logic;
	signal empty_s : std_logic;

	-- write only when there is space, read only when there is data
	signal do_wr : std_logic;
	signal do_rd : std_logic;

	-- registered read data and the depth 1 bypass register
	signal mem_rd_data : std_logic_vector(G_DATA_WIDTH-1 downto 0);
	signal bypass_reg  : std_logic_vector(G_DATA_WIDTH-1 downto 0);

begin

	full_s  <= '1' when count >= i_virtual_depth else '0';
	empty_s <= '1' when count = 0                else '0';

	do_wr <= i_wr_en and not full_s;
	do_rd <= i_rd_en and not empty_s;

	PROC_FIFO_WRITE: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				wr_ptr <= (others => '0');
			elsif do_wr = '1' then
				fifo_mem(to_integer(wr_ptr)) <= i_wr_data;
				if wr_ptr = G_DEPTH-1 then
					wr_ptr <= (others => '0');
				else
					wr_ptr <= wr_ptr + 1;
				end if;
			end if;
		end if;
	end process;

	PROC_FIFO_READ: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				rd_ptr <= (others => '0');
			elsif do_rd = '1' then
				mem_rd_data <= fifo_mem(to_integer(rd_ptr));
				if rd_ptr = G_DEPTH-1 then
					rd_ptr <= (others => '0');
				else
					rd_ptr <= rd_ptr + 1;
				end if;
			end if;
		end if;
	end process;

	-- depth 1 bypass: the memory path needs 2 cycles, so depth 1 uses a plain
	-- register. written on raw i_wr_en, the memory path reads as full at depth 1.
	PROC_BYPASS_REG: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_wr_en = '1' then
				bypass_reg <= i_wr_data;
			end if;
		end if;
	end process;

	o_rd_data <= bypass_reg when i_virtual_depth = 1 else mem_rd_data;

	PROC_FIFO_COUNT: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				count <= (others => '0');
			elsif do_wr = '1' and do_rd = '0' then
				count <= count + 1;
			elsif do_rd = '1' and do_wr = '0' then
				count <= count - 1;
			end if;
		end if;
	end process;

	o_full  <= full_s;
	o_almost_full <= '1' when count >= i_virtual_depth - 1 else '0';
	o_empty <= empty_s;

end architecture rtl;
