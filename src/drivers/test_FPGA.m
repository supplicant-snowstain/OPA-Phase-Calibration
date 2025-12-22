%%
%电源部分
FPGA=serial('com1');%选择设备串口
set(FPGA,'BaudRate',115200,'StopBits',1,'DataBits',8);
fopen(FPGA);
sprintf('FPGA COM ports have been opened.')

ch=1; %通道编号
Volt = 9;%电压单位V
VC96_setVone( FPGA,ch,5*rand(1) );
fclose(FPGA);
delete(instrfindall);
%%
%采集系统部分
FPGA=serial('com1');%选择设备串口
set(FPGA,'BaudRate',115200,'StopBits',1,'DataBits',8);
fopen(FPGA);
sprintf('FPGA COM ports have been opened.')
ch=2;
OS=6;
read_data_pre( FPGA,OS );

for index=1:100
wave(index)=read_data( FPGA,ch );
end
figure;
plot(1:100,wave);
fclose(FPGA);
delete(instrfindall);
