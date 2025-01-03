module accumulatorDCO (
    output [127:0] sum,
    input [127:0] in,
    input clk, 
    input reset, 
    input enable
);
reg [127:0] count;//fo some reason here is 7:0
assign sum=count;

//    always @(posedge clk, posedge reset) begin
//        if (reset) begin
always @ ( posedge clk or negedge reset) begin
	if (~reset) begin
            count <= 0;
        end else if (enable) begin
            count <= count + in;
        end
    end
	

endmodule



