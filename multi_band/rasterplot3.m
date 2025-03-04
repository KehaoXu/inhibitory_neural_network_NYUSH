% % Filename: rasterplot (Version-3)
% % Date: 2025.3.3
% % Author: Jiatong Guo
% % Description: A function to show scatterplots

function [] = rasterplot3(res, param)
ne = param.ne;
ni = param.ni;
ns = param.ns;
nv = param.nv;  % Add VIP neuron count

spike = res.spike;
num_E = sum(spike(1,1:ne));
num_I = sum(spike(1,ne+1:ne+ni));
num_S = sum(spike(1,ne+ni+1:ne+ni+ns));
num_V = sum(spike(1,ne+ni+ns+1:ne+ni+ns+nv));  % Count VIP spikes

coor_E = zeros(num_E,2);
coor_I = zeros(num_I,2);
coor_S = zeros(num_S,2);
coor_V = zeros(num_V,2);  % Coordinates for VIP neurons

index_E = 1;
index_I = 1;
index_S = 1;
index_V = 1;  % Index counter for VIP neurons

ave_E = num_E / ne / param.duration * 1000;
ave_I = num_I / ni / param.duration * 1000;
ave_S = num_S / ns / param.duration * 1000;
ave_V = num_V / nv / param.duration * 1000;  % Average VIP firing rate

% Extract spike coordinates for each neuron type
for i=1:ne
    num_Ei = spike(1,i);
    coor_E(index_E:index_E+num_Ei-1,1) = i;
    coor_E(index_E:index_E+num_Ei-1,2) = spike(2:1+num_Ei,i)*1000;
    index_E = index_E + num_Ei;
end

for i=(ne+1):(ne+ni)
    num_Ii = spike(1,i);
    coor_I(index_I:index_I+num_Ii-1,1) = i;
    coor_I(index_I:index_I+num_Ii-1,2) = spike(2:1+num_Ii,i)*1000;
    index_I = index_I + num_Ii;
end

for i=(ne+ni+1):(ne+ni+ns)
    num_Si = spike(1,i);
    coor_S(index_S:index_S+num_Si-1, 1) = i;
    coor_S(index_S:index_S+num_Si-1, 2) = spike(2:1+num_Si,i)*1000;
    index_S = index_S + num_Si;
end

% Extract VIP spike coordinates
for i=(ne+ni+ns+1):(ne+ni+ns+nv)
    num_Vi = spike(1,i);
    coor_V(index_V:index_V+num_Vi-1, 1) = i;
    coor_V(index_V:index_V+num_Vi-1, 2) = spike(2:1+num_Vi,i)*1000;
    index_V = index_V + num_Vi;
end

% Plot all neuron types with different colors
hE = scatter(coor_E(:,2), coor_E(:,1),10,'.','r');
hold on
hI = scatter(coor_I(:,2), coor_I(:,1),10,'.','b');
hold on
hS = scatter(coor_S(:,2), coor_S(:,1),10,'.','g');
hold on
hV = scatter(coor_V(:,2), coor_V(:,1),10,'.','m');  % VIP neurons in magenta


% Legend & Firing Rate
legend([hE, hI, hS, hV], { ...
    sprintf('PCs (Red): %.2f Hz', ave_E), ...
    sprintf('PVs (Blue): %.2f Hz', ave_I), ...
    sprintf('SOMs (Green): %.2f Hz', ave_S), ...
	sprintf('VIPs (Megenta): %.2f Hz', ave_V)},...
    'FontSize', 15, 'Location', 'southeast');


% Set title and labels
title(sprintf('S_{es} = %.4f, S_{is} = %.4f, \\tau_{es}^{delay} = %d ms, \\tau_{is}^{delay} = %d ms', param.s_es / 100, param.s_is / 100, param.s2e_delay, param.s2i_delay), 'FontSize', 25);
xlabel('time(ms)');
ylabel('Neuron Index');
set(gca,'fontsize',15);
% set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

end