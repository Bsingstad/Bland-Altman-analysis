clear;
d=load('c:\Prosjekter\div\LiDCO vs NexFin\Bland altman kode\sim.txt');
d_zou=struct2array(load('c:\Prosjekter\div\LiDCO vs NexFin\Bland altman kode\zou_tabell.mat'));
%d=d_zou;
s=d(:,1);
ma=d(:,2); mb=d(:,3);
for i=1:length(ma)
    avg_ab(i)=(ma(i)+mb(i))/2;
    diff_ab(i)=ma(i)-mb(i);
end
bias=mean(diff_ab);
std_ab=std(diff_ab);
upper_loa=bias+1.96*sqrt(2*std_ab);
upper_loa=bias+2*std_ab;
lower_loa=bias-2*std_ab;


%Kalkulere mean difference og varians av difference innen hvert individ
subjects=unique(s); n=length(subjects);
for subject=1:n
    idx=find(s==subject);
    mi(subject)=length(idx); %antall målepunkter innen hvert individ
    dmean(subject)=mean(diff_ab(idx));
    dvar(subject)=var(diff_ab(idx));
end

%Overall estimate of random error within subject:
N=sum(mi);
var_dw=sum((mi-1)./(N-n).*dvar);

bias=mean(dmean);

%Estimate between subjects variability:
var_d=var(dmean);

%Harmonic mean av replikater
m_h=n./sum(1./mi);

upper_loa=bias+1.96*sqrt(var_d+(1-1/m_h)*var_dw);
lower_loa=bias-1.96*sqrt(var_d+(1-1/m_h)*var_dw);



%Kalkulere varians for LOA basert på Delta metoden
var_tot=var_d+(1-1/m_h)*var_dw;
var_db=var_tot-var_dw;
var_loa=var_d/n+1.96^2/(2*var_tot)*(var_d^2/(n-1)+(1-1/m_h)^2*(var_dw^2/(N-n)));

upper_loa_upper_ci=upper_loa+1.96*sqrt(var_loa);
upper_loa_lower_ci=upper_loa-1.96*sqrt(var_loa);
lower_loa_upper_ci=lower_loa+1.96*sqrt(var_loa);
lower_loa_lower_ci=lower_loa-1.96*sqrt(var_loa);

bias_upper_ci=bias+tinv(0.975,19)*sqrt(var_d/n); %Finn t alpha/2 fra t tabell
bias_lower_ci=bias-tinv(0.975,19)*sqrt(var_d/n);
%bias_upper_ci=bias+1.96*sqrt(var_tot/n);

figure(1)
gscatter(avg_ab,diff_ab,s)
hold on
plot(avg_ab,repmat(bias,1,length(avg_ab)),'black')
plot(avg_ab,repmat(upper_loa,1,length(avg_ab)),'black')
plot(avg_ab,repmat(lower_loa,1,length(avg_ab)),'black')
plot(avg_ab,repmat(upper_loa_upper_ci,1,length(avg_ab)),'black:')
plot(avg_ab,repmat(upper_loa_lower_ci,1,length(avg_ab)),'black:')
plot(avg_ab,repmat(lower_loa_upper_ci,1,length(avg_ab)),'black:')
plot(avg_ab,repmat(lower_loa_lower_ci,1,length(avg_ab)),'black:')
hold off

%MOVER
l=var_tot-sqrt(((var_d)*(1-(n-1)/(chi2inv(0.975,n-1))))^2 + ((1-1/m_h)*var_dw*(1-(N-n)/(chi2inv(0.975,N-n))))^2);
u=var_tot+sqrt(((var_d)*((n-1)/(chi2inv(0.025,n-1))-1))^2 + ((1-1/m_h)*var_dw*((N-n)/(chi2inv(0.025,N-n))-1))^2);

LME=sqrt(norminv(0.975)^2*(var_d/n)+norminv(0.975)^2*(sqrt(u)-sqrt(var_tot))^2);
RME=sqrt(norminv(0.975)^2*(var_d/n)+norminv(0.975)^2*(sqrt(var_tot)-sqrt(l))^2);

lower_loa_lower_ci_MOVER=lower_loa-LME;
lower_loa_upper_ci_MOVER=lower_loa+RME;
upper_loa_lower_ci_MOVER=upper_loa-RME;
upper_loa_upper_ci_MOVER=upper_loa+LME;