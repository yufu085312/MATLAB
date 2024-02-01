n=10000;
zn=0;
for k=1:n
    zn=(1/(k*k))+zn;
end
%πの近似式
pe=sqrt(6*zn);
%絶対誤差
aerror=abs(pe-pi);
%相対誤差
rerror=abs(aerror/pi);
%10進有効桁数
ddigit=-log10(rerror);

disp(zn)
disp(pe)
disp(aerror)
disp(rerror)
disp(ddigit)