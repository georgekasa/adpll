/*

[Transfer function block diagram]

x[z] --------->| a |------------------->| + |-----> y[z]
         |                                ^
         |                                |
          ----->| b |-----| + |-| z^-1 |---
                            ^             |
                            |             |
                            |             |
                             -------<-----
page 134 book we need an integral loop gain with H(z) = rho * (z^-1)/ (1- z^-1)
*/


`timescale 1ps / 1fs
module theFilter (
    input                         i_clk,                                          // Input clock
    input                         i_rst_an,                                       // Asynchronous active low reset
    input                         i_en,                                           // Registers enable
    input    [3:0]                i_mul_a,                                        // Input multiplier a
    input    [3:0]                i_mul_b,                                        // Input multiplier b
    input   signed [31:0]         i_phase,                                        // Input phase difference
    output  signed [31:0]         o_phase                                         // Output phase
);


   wire signed [31:0] s_mul_a, s_mul_b;
    reg signed [31:0] s_sub_sum_b;
	wire signed [31:0] s_sub_sum_a ;
	assign s_mul_a = i_phase >>> i_mul_a; 
	assign s_mul_b = i_phase >>> i_mul_b;
    assign o_phase = $signed(s_sub_sum_a) + $signed(s_sub_sum_b);

	assign s_sub_sum_a = i_phase >>> i_mul_a; 
    always @(posedge i_clk or negedge i_rst_an) begin
        if(~i_rst_an) begin
            s_sub_sum_b <= 32'd0;
        end
        else if (i_en) begin
            s_sub_sum_b <= $signed(s_sub_sum_b) + s_mul_b;
        end
    end
endmodule