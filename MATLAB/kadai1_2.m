%数値積分
m=[10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000];
Iext=zeros(size(m))+(1/2+log(2));% 厳密解の配列
Imid=zeros(size(m));%中点則
Itrp=zeros(size(m));%台形則
Ismp=zeros(size(m));%シンプソン則
n=size(m,2);
% 分点数を変えて積分
for k=1:n
    Imid(k)=midpoint(@integrand,1,2,m(k)); 
    Itrp(k)=trapezoid(@integrand,1,2,m(k));
    Ismp(k)=simpson(@integrand,1,2,m(k));
end
format shortE;
%　絶対誤差
aemid=abs(Iext-Imid);%中点則
aetrp=abs(Iext-Itrp);%台形則
aesmp=abs(Iext-Ismp);%シンプソン則
% 相対誤差
remid=abs(aemid./Iext);%中点則
retrp=abs(aetrp./Iext);%台形則
resmp=abs(aesmp./Iext);%シンプソン則
%10進有効桁数
ddmidplog=-log10(remid);%中点則
ddtrplog=-log10(retrp);%台形則
ddsmplog=-log10(resmp);%シンプソン則

% 両対数グラフの指定
plot(m, remid, 'r-s'); % 中点則は赤色(r)
hold on;
plot(m, retrp, 'g-s'); % 台形則は緑色(b)
plot(m, resmp, 'b-s'); % シンプソン則は青色(b)
graph=gca;             % グラフハンドルの取得
graph.YScale='log';    % y軸対数スケールの設定
graph.XScale='log';    % x軸対数スケールの設定
xlabel('分点数');
ylabel('相対誤差');
legend('中点則','台形則','シンプソン則');
title('数値積分の誤差')
hold off;

%相対誤差を 10^(-8) 以下とするために必要な分点数
%中点則
mm=(1/2+log(2));
mid=1;
i=1;
while mid > 10^(-8)
   pe=trapezoid(@integrand,1,2,i);
   format short
   mid=abs(pe-mm)/mm;
   i=i+1;
end
disp(i)
%台形則
tm=(1/2+log(2));
trp=1;
h=1;
while trp > 10^(-8)
   pe=midpoint(@integrand,1,2,h);
   format short
   trp=abs(pe-tm)/tm;
   h=h+1;
end
disp(h)
%シンプソン則
sm=(1/2+log(2));
smp=1;
n=1;
while smp > 10^(-8)
   pe=simpson(@integrand,1,2,n);
   format short
   smp=abs(pe-sm)/sm;
   n=n+1;
end
disp(n)

function fout=integrand(x)
    fout=(1+x)/(x*x);
end
% 以下，添付不要
% 中点則　trapezoid
% 関数f(x), 区間[a, b], m個分割
function Mm = midpoint(f,a,b,m)
    h = (b - a) / m;
    ci = a + h/2;
    Mm=0;
    for i = 0:m-1
        Mm = Mm + f(ci);
        ci = ci + h;
    end
    Mm=Mm*h;
end

% 台形則　trapezoid
% 関数f(x), 区間[a, b], m個分割
function Tm = trapezoid (f,a,b,m)
    h = (b - a) / m;
    ai = a + h;
    Tm=(f(a)+f(b)) / 2;
    for i = 1:m-1
        Tm = Tm + f(ai);
        ai = ai + h;
    end
    Tm=Tm*h;
end

% シンプソン則　simpson
% 関数f(x), 区間[a, b], m個分割
function Sm = simpson (f,a,b,m)
    % 中点則と台形則を計算
    Mm= midpoint(f, a, b, m/2);
    Tm= trapezoid(f, a, b, m/2);
    % 中点則と台形則の結果を合成
    Sm = (2*Mm + Tm) / 3;
end