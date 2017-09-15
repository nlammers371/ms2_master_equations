%set size of PolII Blocking distance in bp
b_s = [50,100,150];
%define vector of elongation rates (bp/s)
v_e = (1:.2:4).*1000/60;
%define range of arrival rates to test
r_a = linspace(.01,5);
%set reference threshold level
p = .1;
t_fig = figure;
colormap('winter');
cm = colormap;
increment = floor(60/(1+length(v_e)));
for i = 1:length(b_s)
    b = b_s(i);
    subplot(length(b_s),1,i); 
    hold on
    for j = 1:length(v_e)
        v = v_e(j);
        factor = sqrt(b/v^2 + 1./r_a.^2)./ (b/v + 1./r_a);
        plot(r_a,factor,'Color', cm(increment*j,:))
    end
    plot(r_a,linspace(b^-.5,b^-.5),'--','Color',cm(1,:))
    title(['Coefficient of Variation: Time per Loading Event (PolII Size: ' num2str(b) ')']);
    hold off
end
saveas(t_fig,'../figs/t_coeff_var.png','png')