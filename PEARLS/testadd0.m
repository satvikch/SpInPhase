n=0:127;

f=0.1;

y1=sin(2*pi*f*n);

y1(1:16)=y1(1:16)+1;


y2s=y1(1:32);

fy2s=abs(fft(y2s));

y2=[y2s zeros(1,128-32)];

y3=[zeros(1,128-32),y2s];

figure;

plot(abs(fft(y1)));


hold on;


plot(abs(fft(y2)),'r');
plot(abs(fft(y3)),'g');
stem(1:128/32:128,fy2s,'k');
hold off;



figure;

plot(angle(fft(y1)));
hold on;
plot(angle(fft(y2)),'r');
plot(angle(fft(y3)),'g');

hold off;