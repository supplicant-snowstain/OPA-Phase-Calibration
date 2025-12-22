function Volt = VC96_setVone( spd,chv,Volt )
%第8位不取反，完全按照协议
% CH_num = 1;% ch编号与PCB板标号之间的对应关系见: DB25引脚配置.xlsx
% Volt = 9;%电压单位V
[DAC_num, CH_num]=VC96_getnum(chv);
strinit(1) = hex2dec('ff');
strinit(2) = hex2dec('00')+DAC_num;
strinit(3) = hex2dec('7b');
strinit(4) = 0;
strinit(5) = 1;
strinit(6) = rem((strinit(1)+strinit(2)+strinit(3)+strinit(4)+strinit(5)),256);
fwrite(spd,strinit);

strinit(1) = hex2dec('ff');
strinit(2) = hex2dec('00')+DAC_num;
strinit(3) = CH_num + (hex2dec('04')*16);
temp=floor(Volt/10*65536);
string=dec2bin(temp,16);
strinit(4) =bin2dec(string(1:8));
strinit(5) =bin2dec(string(9:16));
strinit(6) = rem((strinit(1)+strinit(2)+strinit(3)+strinit(4)+strinit(5)),256);
tic
fwrite(spd,strinit);
toc
end

function [DAC_num, CH_num]=VC96_getnum(chv)
%For ease of understanding, the parameter chv starts from 1
chv=chv-1; %real chanel number starts from 0, but not 1.
DAC_num_order=floor(chv/16)+1;
DAC_num_list=[0, 8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15];
DAC_num=DAC_num_list(DAC_num_order);
CH_num=mod(chv,16);
end

