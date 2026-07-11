library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

-- Self-checking testbench for the synchronous FIFO.
-- Exercises: reset, fill-to-full, overflow protection, drain-to-empty,
-- underflow protection, read-data ordering/latency, and pointer wrap-around.
-- Tool-agnostic VHDL-2008 (assert/report/severity); runs under GHDL and xsim.
entity tb_fifo is
end entity tb_fifo;

architecture sim of tb_fifo is

    constant C_DW    : integer := 8;   -- data width
    constant C_DEPTH : integer := 5;   -- FIFO depth (non-power-of-two on purpose)
    constant C_PERIOD : time := 10 ns;

    signal clk     : std_logic := '0';
    signal reset   : std_logic := '1';
    signal wr_en   : std_logic := '0';
    signal wr_data : std_logic_vector(C_DW-1 downto 0) := (others => '0');
    signal rd_en   : std_logic := '0';
    signal rd_data : std_logic_vector(C_DW-1 downto 0);
    signal full    : std_logic;
    signal empty   : std_logic;

    signal done    : boolean := false;
    signal errors  : integer := 0;

    -- Compare a std_logic_vector against an expected integer and count mismatches.
    procedure check_data(signal err : inout integer;
                         got : std_logic_vector; exp : integer; msg : string) is
    begin
        if to_integer(unsigned(got)) /= exp then
            report "FAIL: " & msg & " got=" & integer'image(to_integer(unsigned(got))) &
                   " exp=" & integer'image(exp) severity error;
            err <= err + 1;
        end if;
    end procedure;

    procedure check_flag(signal err : inout integer;
                         got : std_logic; exp : std_logic; msg : string) is
    begin
        if got /= exp then
            report "FAIL: " & msg & " got=" & std_logic'image(got) &
                   " exp=" & std_logic'image(exp) severity error;
            err <= err + 1;
        end if;
    end procedure;

begin

    dut: entity work.fifo
        generic map (
            g_data_width => C_DW,
            g_depth      => C_DEPTH
        )
        port map (
            i_clk     => clk,
            i_reset   => reset,
            i_wr_en   => wr_en,
            i_wr_data => wr_data,
            i_rd_en   => rd_en,
            o_rd_data => rd_data,
            o_full    => full,
            o_empty   => empty
        );

    clk <= not clk after C_PERIOD/2 when not done else '0';

    stimulus: process
        -- Sample a combinational/registered value shortly after the clock edge.
        procedure settle is
        begin
            wait for C_PERIOD/10;
        end procedure;
    begin
        --------------------------------------------------------------------
        -- 1. Reset behaviour
        --------------------------------------------------------------------
        reset <= '1';
        wait until rising_edge(clk);
        wait until rising_edge(clk);
        reset <= '0';
        settle;
        check_flag(errors, empty, '1', "after reset FIFO should be empty");
        check_flag(errors, full,  '0', "after reset FIFO should not be full");

        --------------------------------------------------------------------
        -- 2. Fill the FIFO with C_DEPTH known words (10, 20, 30, 40)
        --------------------------------------------------------------------
        for i in 0 to C_DEPTH-1 loop
            wr_en   <= '1';
            wr_data <= std_logic_vector(to_unsigned((i+1)*10, C_DW));
            wait until rising_edge(clk);
        end loop;
        wr_en <= '0';
        settle;
        check_flag(errors, full,  '1', "FIFO should be full after C_DEPTH writes");
        check_flag(errors, empty, '0', "FIFO should not be empty when full");

        --------------------------------------------------------------------
        -- 3. Overflow protection: extra write while full must be ignored
        --------------------------------------------------------------------
        wr_en   <= '1';
        wr_data <= std_logic_vector(to_unsigned(99, C_DW));
        wait until rising_edge(clk);
        wr_en <= '0';
        settle;
        check_flag(errors, full, '1', "FIFO must stay full after blocked overflow write");

        --------------------------------------------------------------------
        -- 4. Drain the FIFO and check ordering + 1-cycle read latency.
        --    o_rd_data becomes valid on the edge that consumes the word.
        --------------------------------------------------------------------
        for i in 0 to C_DEPTH-1 loop
            rd_en <= '1';
            wait until rising_edge(clk);   -- edge consumes word i, registers rd_data
            rd_en <= '0';
            settle;
            check_data(errors, rd_data, (i+1)*10,
                       "read data mismatch at index " & integer'image(i));
            wait until rising_edge(clk);   -- gap cycle between single reads
        end loop;
        settle;
        check_flag(errors, empty, '1', "FIFO should be empty after draining");
        check_flag(errors, full,  '0', "FIFO should not be full when empty");

        --------------------------------------------------------------------
        -- 5. Underflow protection: read while empty must not corrupt state
        --------------------------------------------------------------------
        rd_en <= '1';
        wait until rising_edge(clk);
        rd_en <= '0';
        settle;
        check_flag(errors, empty, '1', "FIFO must stay empty after blocked underflow read");

        --------------------------------------------------------------------
        -- 6. Wrap-around: push past the address-wrap boundary and drain.
        --    Write 3, read 3, then write 4 (wrapping the address) and read 4.
        --------------------------------------------------------------------
        for i in 0 to 2 loop
            wr_en   <= '1';
            wr_data <= std_logic_vector(to_unsigned(100 + i, C_DW));
            wait until rising_edge(clk);
        end loop;
        wr_en <= '0';

        for i in 0 to 2 loop
            rd_en <= '1';
            wait until rising_edge(clk);
            rd_en <= '0';
            settle;
            check_data(errors, rd_data, 100 + i,
                       "wrap phase-1 read mismatch at index " & integer'image(i));
            wait until rising_edge(clk);
        end loop;

        for i in 0 to C_DEPTH-1 loop
            wr_en   <= '1';
            wr_data <= std_logic_vector(to_unsigned(200 + i, C_DW));
            wait until rising_edge(clk);
        end loop;
        wr_en <= '0';
        settle;
        check_flag(errors, full, '1', "FIFO should be full after wrap-around fill");

        for i in 0 to C_DEPTH-1 loop
            rd_en <= '1';
            wait until rising_edge(clk);
            rd_en <= '0';
            settle;
            check_data(errors, rd_data, 200 + i,
                       "wrap phase-2 read mismatch at index " & integer'image(i));
            wait until rising_edge(clk);
        end loop;
        settle;
        check_flag(errors, empty, '1', "FIFO should be empty after wrap-around drain");

        --------------------------------------------------------------------
        -- Report result
        --------------------------------------------------------------------
        if errors = 0 then
            report "TB PASS -- all FIFO checks passed" severity note;
        else
            report "TB FAIL -- " & integer'image(errors) & " check(s) failed" severity failure;
        end if;

        done <= true;
        wait;
    end process stimulus;

end architecture sim;
