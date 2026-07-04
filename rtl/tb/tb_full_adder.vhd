library ieee;
use ieee.std_logic_1164.all;

entity tb_full_adder is
end entity tb_full_adder;

architecture sim of tb_full_adder is
    signal a    : std_logic := '0';
    signal b    : std_logic := '0';
    signal cin  : std_logic := '0';
    signal sum  : std_logic;
    signal cout : std_logic;
begin
    dut: entity work.full_adder
        port map (
            a    => a,
            b    => b,
            cin  => cin,
            sum  => sum,
            cout => cout
        );

    stimulus: process
    begin
        a <= '0'; b <= '0'; cin <= '0'; wait for 10 ns;
        a <= '0'; b <= '0'; cin <= '1'; wait for 10 ns;
        a <= '0'; b <= '1'; cin <= '0'; wait for 10 ns;
        a <= '0'; b <= '1'; cin <= '1'; wait for 10 ns;
        a <= '1'; b <= '0'; cin <= '0'; wait for 10 ns;
        a <= '1'; b <= '0'; cin <= '1'; wait for 10 ns;
        a <= '1'; b <= '1'; cin <= '0'; wait for 10 ns;
        a <= '1'; b <= '1'; cin <= '1'; wait for 10 ns;

        wait;
    end process stimulus;
end architecture sim;