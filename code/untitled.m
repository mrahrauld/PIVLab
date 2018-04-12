close all
[D , V]  = textread('line_genval.txt','%f %f','delimiter',',');
num = xlsread('true_genval.xlsx');
Title = 'Vitesses de surface de l''Argentine';
plot_graph(D,V,num(:,1),num(:,2),Title);


[D , V]  = textread('line_wimbe.txt','%f %f','delimiter',',');
num = xlsread('true_wimbe.xlsx');
Title = 'Vitesses de surface de la Wimbe';
plot_graph(D,V,num(:,1),num(:,2),Title);

[D , V]  = textread('line_daverdisse1.txt','%f %f','delimiter',',');
num = xlsread('true_daverdisse1.xlsx');
Title = 'Vitesses de surface, Daverdisse';
plot_graph(D,V,num(:,1),num(:,2),Title);

[D , V]  = textread('line_eprave_pont.txt','%f %f','delimiter',',');
num = xlsread('true_eprave.xlsx');
Title = 'Vitesses de surface de l''Eprave';
plot_graph(D,V,num(:,1),num(:,2),Title);

[D , V]  = textread('line_eprave_eau.txt','%f %f','delimiter',',');
num = xlsread('true_eprave.xlsx');
Title = 'Vitesses de surface de l''Eprave';
plot_graph(D,V,num(:,1),num(:,2),Title);

function [] = plot_graph(D_LSPIV,V_LSPIV,D_TRUE,V_TRUE,Title)
    D_LSPIV = D_LSPIV-D_LSPIV(1);
    D_TRUE = D_TRUE-D_TRUE(1);
    D_TRUE = D_TRUE/D_TRUE(end);
    D_TRUE = D_TRUE*D_LSPIV(end);
    figure;
    linewidth = 2;
    p = plot(D_LSPIV,V_LSPIV,'r',D_TRUE,V_TRUE,'b--');
    p(1).LineWidth = linewidth;
    p(2).LineWidth = linewidth;
    title(Title,'fontsize',18);
    xlabel('distance (m)','fontsize',14) % x-axis label
    ylabel('vitesse (m/s)','fontsize',14) % y-axis label
    lgd = legend('LSPIV','ADC','Location','northwest');
    lgd.FontSize = 14;
    MSE = immse(V_LSPIV,spline(D_TRUE,V_TRUE,D_LSPIV));
    MSE2 = immse(V_TRUE,spline(D_LSPIV,V_LSPIV,D_TRUE));
    idx = find(V_TRUE);
    inter = abs(spline(D_LSPIV,V_LSPIV,D_TRUE)-V_TRUE)./V_TRUE;
    Mean_perc_error = mean(inter(idx))*100
    %inter2 = sum(inter(idx).*V_TRUE(idx))/sum(V_TRUE(idx))
end
