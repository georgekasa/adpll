
`timescale 1ps / 1fs
module dff(
in   , // Data Input
out , //output
clk    , // Clock Input
reset  //, // Reset input

 );
//-----------Input Ports---------------
input in, clk, reset ; 

//-----------Output Ports---------------
output out;

//------------Internal Variables--------
reg q;
assign out = q;

//-------------Code Starts Here---------
always @ ( posedge clk or negedge reset) begin
	if (~reset) begin
 		 q <= 1'd0;
	end  else begin
  		q <= in;
	end

end

endmodule