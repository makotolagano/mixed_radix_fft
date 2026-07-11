library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity fifo is
	generic (
		g_data_width : integer := 36;
		g_depth 		 : integer := 1024
	);
	port (
		i_clk 		: in std_logic;
		i_reset 	: in std_logic;

		i_wr_en 	: in std_logic;
		i_wr_data : in std_logic_vector(g_data_width-1 downto 0);

		i_rd_en 	: in std_logic;
		o_rd_data : out std_logic_vector(g_data_width-1 downto 0);

		o_full 		: out std_logic;
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

	constant C_ADDR_W : integer := clogb2(g_depth);      -- address bits
	constant C_CNT_W  : integer := clogb2(g_depth + 1);  -- occupancy 0..g_depth

	-- Type definitions
	type t_mem is array (0 to g_depth-1) of std_logic_vector(g_data_width-1 downto 0);
	signal fifo_mem : t_mem;

	-- Pointers wrap at g_depth (not necessarily a power of two), so an
	-- occupancy counter is used to distinguish full from empty.
	signal wr_ptr : unsigned(C_ADDR_W-1 downto 0) := (others => '0');
	signal rd_ptr : unsigned(C_ADDR_W-1 downto 0) := (others => '0');
	signal count  : unsigned(C_CNT_W-1 downto 0)  := (others => '0');

	signal full_s  : std_logic;
	signal empty_s : std_logic;

	-- Qualified enables (write only when space, read only when data).
	signal do_wr : std_logic;
	signal do_rd : std_logic;

begin

	full_s  <= '1' when count = g_depth else '0';
	empty_s <= '1' when count = 0       else '0';

	do_wr <= i_wr_en and not full_s;
	do_rd <= i_rd_en and not empty_s;

	PROC_FIFO_WRITE: process(i_clk)
	begin
		if rising_edge(i_clk) then
			if i_reset = '1' then
				wr_ptr <= (others => '0');
			elsif do_wr = '1' then
				fifo_mem(to_integer(wr_ptr)) <= i_wr_data;
				if wr_ptr = g_depth-1 then
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
				o_rd_data <= fifo_mem(to_integer(rd_ptr));
				if rd_ptr = g_depth-1 then
					rd_ptr <= (others => '0');
				else
					rd_ptr <= rd_ptr + 1;
				end if;
			end if;
		end if;
	end process;

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
	o_empty <= empty_s;

end architecture rtl;
