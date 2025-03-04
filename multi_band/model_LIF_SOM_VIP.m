% % Filename: model_LIF_SOM_VIP
% % Date: 2025.3.3
% % Author: Jiatong Guo
% % Description: Add VIP neurons


function [res] = model_LIF_SOM_VIP(param, init)

%% Parameter sets 
% e:excitatory, i:inhibitory, s:somatostatin, v:VIP
ne      = param.ne;     % number of excitatory cells
ni      = param.ni;     % number of PV cells
ns      = param.ns;     % number of SOM cells
nv      = param.nv;     % number of VIP cells

p_ee    = param.p_ee;   % projection probability
p_ie    = param.p_ie;
p_ei    = param.p_ei;
p_ii    = param.p_ii;
p_se    = param.p_se;
p_es    = param.p_es;
p_is    = param.p_is;
p_si    = param.p_si;
p_sv    = param.p_sv;   % from SOM to VIP
p_vs    = param.p_vs;   % from VIP to SOM
p_ve    = param.p_ve;   % from VIP to Excitatory
p_vi    = param.p_vi;   % from VIP to PV
p_ev    = param.p_ev;   % from Excitatory to VIP
p_iv    = param.p_iv;   % from PV to VIP
p_vv    = param.p_vv;   % from VIP to VIP (typically low)

s_ee    = param.s_ee;   % synaptic strength
s_ie    = param.s_ie;
s_ei    = param.s_ei;
s_ii    = param.s_ii;
s_se    = param.s_se;
s_is    = param.s_is;
s_si    = param.s_si;
s_es    = param.s_es;
s_sv    = param.s_sv;   % from SOM to VIP
s_vs    = param.s_vs;   % from VIP to SOM
s_ve    = param.s_ve;   % from VIP to Excitatory
s_vi    = param.s_vi;   % from VIP to PV
s_ev    = param.s_ev;   % from Excitatory to VIP
s_iv    = param.s_iv;   % from PV to VIP
s_vv    = param.s_vv;   % from VIP to VIP

s_exe   = param.s_exe;  % external input strength
s_exi   = param.s_exi;
s_exs   = param.s_exs;
s_exv   = param.s_exv;  % external input strength to VIP

tau_ee  = param.tau_ee;     % time constant
tau_ie  = param.tau_ie;
tau_ei  = param.tau_ei;
tau_ii  = param.tau_ii;
tau_se  = param.tau_se;
tau_is  = param.tau_is;
tau_si  = param.tau_si;
tau_es  = param.tau_es;
tau_ev  = param.tau_ev;     % Excitatory to VIP
tau_iv  = param.tau_iv;     % PV to VIP
tau_sv  = param.tau_sv;     % SOM to VIP
tau_ve  = param.tau_ve;     % VIP to Excitatory
tau_vi  = param.tau_vi;     % VIP to PV
tau_vs  = param.tau_vs;     % VIP to SOM
tau_vv  = param.tau_vv;     % VIP to VIP

dt      = param.gridsize;
tau_re  = param.tau_re;     % refractory period for excitatory neurons
tau_ri  = param.tau_ri;     % refractory period for PV
tau_rs  = param.tau_rs;     % refractory period for SOM
tau_rv  = param.tau_rv;     % refractory period for VIP

M       = param.M;      % spike threshold
Mr      = param.Mr;     % inhibitory reversal potential
duration = param.duration; % ms

lambda_e = param.lambda_e/1000; % external spike frequency
lambda_i = param.lambda_i/1000;
lambda_s = param.lambda_s/1000;
lambda_v = param.lambda_v/1000; % external spike frequency for VIP

% SFA parameters for SOM and VIP neurons
tau_as   = param.tau_as;    % adaptation time constant for SOM
b_s      = param.b_s;       % adaptation strength for SOM
tau_av   = param.tau_av;    % adaptation time constant for VIP
b_v      = param.b_v;       % adaptation strength for VIP

delay           = 2;
flag_off_time   = 0.0;
start_threshold = 2;
end_threshold   = 2;
MFE_interval    = 1;
flag_time       = 0;

% time delay of dendritic inhibition
s2e_delay = param.s2e_delay;            % time delay for dendritic inhibition from SOM cells
s2i_delay = param.s2i_delay;            % time delay for inhibition from SOM to PV
s2e_delay_steps = ceil(s2e_delay / dt);    % transfer time delay to steps
s2i_delay_steps = ceil(s2i_delay / dt);
s2e_inhibition_buffer = zeros(s2e_delay_steps, ne);
s2i_inhibition_buffer = zeros(s2i_delay_steps, ni);

% generate connection between neurons
connection_matrix_e = zeros(ne, ne+ni+ns+nv);
connection_matrix_i = zeros(ni, ne+ni+ns+nv);
connection_matrix_s = zeros(ns, ne+ni+ns+nv);
connection_matrix_v = zeros(nv, ne+ni+ns+nv);

% Set connections based on probabilities
connection_matrix_e(:,1:ne)                 = binornd(1,p_ee,ne,ne);      % E to E
connection_matrix_i(:,1:ne)                 = binornd(1,p_ei,ni,ne);      % E to PV
connection_matrix_s(:,1:ne)                 = binornd(1,p_es,ns,ne);      % E to SOM
connection_matrix_v(:,1:ne)                 = binornd(1,p_ev,nv,ne);      % E to VIP

connection_matrix_e(:,ne+1:ne+ni)           = binornd(1,p_ie,ne,ni);      % PV to E
connection_matrix_i(:,ne+1:ne+ni)           = binornd(1,p_ii,ni,ni);      % PV to PV
connection_matrix_s(:,ne+1:ne+ni)           = binornd(1,p_is,ns,ni);      % PV to SOM
connection_matrix_v(:,ne+1:ne+ni)           = binornd(1,p_iv,nv,ni);      % PV to VIP

connection_matrix_e(:,ne+ni+1:ne+ni+ns)     = binornd(1,p_se,ne,ns);      % SOM to E
connection_matrix_i(:,ne+ni+1:ne+ni+ns)     = binornd(1,p_si,ni,ns);      % SOM to PV
connection_matrix_s(:,ne+ni+1:ne+ni+ns)     = binornd(1,0,ns,ns);         % SOM to SOM (0 per Hertäg paper)
connection_matrix_v(:,ne+ni+1:ne+ni+ns)     = binornd(1,p_sv,nv,ns);      % SOM to VIP

connection_matrix_e(:,ne+ni+ns+1:ne+ni+ns+nv) = binornd(1,p_ve,ne,nv);    % VIP to E
connection_matrix_i(:,ne+ni+ns+1:ne+ni+ns+nv) = binornd(1,p_vi,ni,nv);    % VIP to PV
connection_matrix_s(:,ne+ni+ns+1:ne+ni+ns+nv) = binornd(1,p_vs,ns,nv);    % VIP to SOM
connection_matrix_v(:,ne+ni+ns+1:ne+ni+ns+nv) = binornd(1,p_vv,nv,nv);    % VIP to VIP (typically near 0)

% Combine all connection matrices and eliminate self-connections
connection_mat = [connection_matrix_e; connection_matrix_i; connection_matrix_s; connection_matrix_v];
connection_mat(logical(eye(ne+ni+ns+nv))) = 0;

% Extract individual connection matrices
connection_matrix_e = connection_mat(1:ne,:);
connection_matrix_i = connection_mat(ne+1:ne+ni,:);
connection_matrix_s = connection_mat(ne+ni+1:ne+ni+ns,:);
connection_matrix_v = connection_mat(ne+ni+ns+1:ne+ni+ns+nv,:);

% determine External Spike Interval according to lambda
% row: external spike;  column: neuron index
esi_e = exprnd(1/lambda_e,ceil(lambda_e*(duration*1.4)),ne); 
esi_i = exprnd(1/lambda_i,ceil(lambda_i*(duration*1.4)),ni);
esi_s = exprnd(1/lambda_s,ceil(lambda_s*(duration*1.4)),ns);
esi_v = exprnd(1/lambda_v,ceil(lambda_v*(duration*1.4)),nv);

while min(sum(esi_e))<duration || min(sum(esi_i))<duration || min(sum(esi_s))<duration || min(sum(esi_v))<duration
   esi_e = exprnd(1/lambda_e,ceil(lambda_e*(duration*1.4)),ne); 
   esi_i = exprnd(1/lambda_i,ceil(lambda_i*(duration*1.4)),ni); 
   esi_s = exprnd(1/lambda_s,ceil(lambda_s*(duration*1.4)),ns);
   esi_v = exprnd(1/lambda_v,ceil(lambda_v*(duration*1.4)),nv);
end

%% convert esi to the external effect at each time point.
ex_e = zeros(duration/dt,ne);
ex_i = zeros(duration/dt,ni);
ex_s = zeros(duration/dt,ns);
ex_v = zeros(duration/dt,nv);

for i=1:ne      % neuron index
   t=esi_e(1,i);
   count=1;    % spike index
   while t<duration
       ind = ceil(t/dt);
       ex_e(ind,i) = ex_e(ind,i) + exp((t-ind*dt)/tau_ee)/tau_ee;
       count=count+1;
       t = t+esi_e(count,i);
   end
end

for i=1:ni
   t=esi_i(1,i);
   count=1;
   while t<duration
       ind=ceil(t/dt);
       ex_i(ind,i)=ex_i(ind,i)+exp((t-ind*dt)/tau_ie)/tau_ie;
       count=count+1;
       t=t+esi_i(count,i);
   end
end

for i=1:ns
   t=esi_s(1,i);
   count=1;
   while t<duration
       ind=ceil(t/dt);
       ex_s(ind,i)=ex_s(ind,i)+exp((t-ind*dt)/tau_se)/tau_se;
       count=count+1;
       t=t+esi_s(count,i);
   end
end

for i=1:nv
   t=esi_v(1,i);
   count=1;
   while t<duration
       ind=ceil(t/dt);
       ex_v(ind,i)=ex_v(ind,i)+exp((t-ind*dt)/tau_ev)/tau_ev;
       count=count+1;
       t=t+esi_v(count,i);
   end
end

qe=exp(-dt/tau_ee);
qi=exp(-dt/tau_ie);
qs=exp(-dt/tau_se);
qv=exp(-dt/tau_ev);

for i=2:size(ex_e,1)
   ex_e(i,:)=ex_e(i,:)+ex_e(i-1,:)*qe; % consider the smearing effect of each spike
   ex_i(i,:)=ex_i(i,:)+ex_i(i-1,:)*qi;
   ex_s(i,:)=ex_s(i,:)+ex_s(i-1,:)*qs;
   ex_v(i,:)=ex_v(i,:)+ex_v(i-1,:)*qv;
end

if isempty(init)
   ve=zeros(1,ne); % membrane potential
   vi=zeros(1,ni);
   vs=zeros(1,ns);
   vv=zeros(1,nv);
   he=zeros(3,ne); % row1: excitation; row2: PV inhibition; row3: SOM inhibition (dendritic)
   hi=zeros(3,ni); % row1: excitation; row2: PV inhibition; row3: SOM inhibition
   hs=zeros(3,ns); % row1: excitation; row2: PV inhibition; row3: VIP inhibition
   hv=zeros(3,nv); % row1: excitation; row2: PV inhibition; row3: SOM inhibition
   % Adaptation variables
   as=zeros(1,ns); % adaptation for SOM
   av=zeros(1,nv); % adaptation for VIP
else
   ve=init.ve;
   vi=init.vi;
   vs=init.vs;
   vv=init.vv;
   he=init.he;
   hi=init.hi;
   hs=init.hs;
   hv=init.hv;
   as=init.as;
   av=init.av;
end

res.HE          = zeros(ceil(duration/dt)+1,ne+ni+ns+nv);
res.HI          = zeros(ceil(duration/dt)+1,ne+ni+ns+nv);
res.HS          = zeros(ceil(duration/dt)+1,ne+ni+ns+nv);
res.VE          = zeros(ceil(duration/dt)+1,ne);
res.VI          = zeros(ceil(duration/dt)+1,ni);
res.VS          = zeros(ceil(duration/dt)+1,ns);
res.VV          = zeros(ceil(duration/dt)+1,nv);
res.AS          = zeros(ceil(duration/dt)+1,ns); % Adaptation for SOM
res.AV          = zeros(ceil(duration/dt)+1,nv); % Adaptation for VIP
res.spikecount_e= zeros(ceil(duration/dt)+1,1);
res.spikecount_i= zeros(ceil(duration/dt)+1,1);
res.spikecount_s= zeros(ceil(duration/dt)+1,1);
res.spikecount_v= zeros(ceil(duration/dt)+1,1);
res.MFE_time    = zeros(ceil(duration)+1,   2);
res.HEE_stat    = zeros(ceil(duration)+1,   3);

MFE_time2       = zeros(ceil(duration)+1,2);
HEE_stat2       = zeros(ceil(duration)+1,3);
wave_spike_count2 = zeros(1,ceil(duration)+1);
wave_spike_count = zeros(1,ceil(duration)+1);
wave_record     = zeros(1,100000);
wave_count      = 1;

refc_e = zeros(1,ne);   % reference clock, represents the time since last spike
refc_i = zeros(1,ni);
refc_s = zeros(1,ns);
refc_v = zeros(1,nv);
sl     = duration;

spike_row = sl * 10;
spike_e = zeros(spike_row,ne);
spike_i = zeros(spike_row,ni);
spike_s = zeros(spike_row,ns);
spike_v = zeros(spike_row,nv);
nearray = 0:ne-1;
niarray = 0:ni-1;
nsarray = 0:ns-1;
nvarray = 0:nv-1;

flag = 0;   % wave data processing flag

for step = 2:duration/dt
   time = step*dt;

   % Handle refractory periods
   rind_e = abs(refc_e)<10^-7; % reference index, if true, the neuron can spike, else in refractory period
   rind_i = abs(refc_i)<10^-7;
   rind_s = abs(refc_s)<10^-7;
   rind_v = abs(refc_v)<10^-7;
   
   refc_e(~rind_e) = refc_e(~rind_e) - dt; % for neurons still in refractory period
   refc_i(~rind_i) = refc_i(~rind_i) - dt;
   refc_s(~rind_s) = refc_s(~rind_s) - dt;
   refc_v(~rind_v) = refc_v(~rind_v) - dt;
   
   % Membrane potential update with adaptation
   ve(rind_e) = ve(rind_e) + dt * (ex_e(step,rind_e)*s_exe  + s_ee*he(1,rind_e)/tau_ee - s_ei*he(2,rind_e)/tau_ei.*(ve(rind_e)+Mr)/(M+Mr)	- s_es*he(3,rind_e)/tau_es.*(ve(rind_e)+Mr)/(M+Mr));
   vi(rind_i)	= vi(rind_i)	+ dt * (ex_i(step,rind_i)*s_exi		+ s_ie*hi(1,rind_i)/tau_ie		- s_ii*hi(2,rind_i)/tau_ii.*(vi(rind_i)+Mr)/(M+Mr)			- s_is*hi(3,rind_i)/tau_is.*(vi(rind_i)+Mr)/(M+Mr));
   vs(rind_s) = vs(rind_s)	+ dt * (ex_s(step,rind_s)*s_exs	+ s_se*hs(1,rind_s)/tau_se - s_si*hs(2,rind_s)/tau_si.*(vs(rind_s)+Mr)/(M+Mr)		- s_vs*hs(3,rind_s)/tau_vs.*(vs(rind_s)+Mr)/(M+Mr) - as(rind_s));
   vv(rind_v) = vv(rind_v)	+ dt * (ex_v(step,rind_v)*s_exv   + s_ev*hv(1,rind_v)/tau_ev - s_iv*hv(2,rind_v)/tau_iv.*(vv(rind_v)+Mr)/(M+Mr)		- s_sv*hv(3,rind_v)/tau_sv.*(vv(rind_v)+Mr)/(M+Mr) - av(rind_v));
   
   % Exponential decay of input in each step
   he = he.*[exp(-dt/tau_ee); exp(-dt/tau_ei); exp(-dt/tau_es)];   
   hi = hi.*[exp(-dt/tau_ie); exp(-dt/tau_ii); exp(-dt/tau_is)];
   hs = hs.*[exp(-dt/tau_se); exp(-dt/tau_si); exp(-dt/tau_vs)];
   hv = hv.*[exp(-dt/tau_ev); exp(-dt/tau_iv); exp(-dt/tau_sv)];
   
   % Update adaptation variables
   as = as + dt * (-as/tau_as + b_s * (vs > M));
   av = av + dt * (-av/tau_av + b_v * (vv > M));
   
   % Spike detection
   sind_e = ve > M; 
   sind_i = vi > M; 
   sind_s = vs > M;
   sind_v = vv > M;
   
   spikecount_e = sum(sind_e);
   spikecount_i = sum(sind_i);
   spikecount_s = sum(sind_s);
   spikecount_v = sum(sind_v);

   % Spike tracking for wave detection
   if flag==1
       wave_spike_count(wave_count) = wave_spike_count(wave_count) + spikecount_e + spikecount_i + spikecount_s + spikecount_v; 
   end
   
   wave_record(wave_record(1)+2 : wave_record(1)+ spikecount_e+1) = time;
   wave_record(1) = wave_record(1) + spikecount_e;
   
   % Process excitatory spikes
   if spikecount_e > 0
       if tau_re < 0  
           ve(sind_e) = ve(sind_e) - M;    
       else
           ve(sind_e) = 0; % reset to zero
           refc_e(sind_e) = refc_e(sind_e) + tau_re;    % update reference clock
       end

       % Update post-synaptic inputs for excitatory spikes
       he(1,:) = he(1,:) + sum(connection_matrix_e(sind_e, 1:ne), 1);    
       hi(1,:) = hi(1,:) + sum(connection_matrix_e(sind_e, ne+1:ne+ni), 1);
       hs(1,:) = hs(1,:) + sum(connection_matrix_e(sind_e, ne+ni+1:ne+ni+ns), 1);
       hv(1,:) = hv(1,:) + sum(connection_matrix_e(sind_e, ne+ni+ns+1:ne+ni+ns+nv), 1);

       % Record spike times
       spike_e(1,sind_e) = spike_e(1,sind_e) + 1;  % spike number + 1 for each firing neuron              
       spikeind = nearray(sind_e)*spike_row + spike_e(1,sind_e) + 1;
       spike_e(round(spikeind)) = step*dt;
   end

   % Process PV (inhibitory) spikes
   if spikecount_i > 0
        if tau_ri < 0
           vi(sind_i) = vi(sind_i) - M;
        else
           vi(sind_i) = 0;
           refc_i(sind_i) = refc_i(sind_i) + tau_ri;
        end
       
       % Update post-synaptic inputs for PV spikes
       he(2,:) = he(2,:) + sum(connection_matrix_i(sind_i, 1:ne), 1);
       hi(2,:) = hi(2,:) + sum(connection_matrix_i(sind_i, ne+1:ne+ni), 1);
       hs(2,:) = hs(2,:) + sum(connection_matrix_i(sind_i, ne+ni+1:ne+ni+ns), 1);
       hv(2,:) = hv(2,:) + sum(connection_matrix_i(sind_i, ne+ni+ns+1:ne+ni+ns+nv), 1);

       % Record spike times
       spike_i(1,sind_i) = spike_i(1,sind_i) + 1; 
       spikeind = niarray(sind_i)*spike_row + spike_i(1,sind_i) + 1;
       spike_i(round(spikeind)) = step*dt;
   end

   % Obtain delayed inhibition from SOM to PC and PV
   he(3,:) = he(3,:) + s2e_inhibition_buffer(1, :);
   hi(3,:) = hi(3,:) + s2i_inhibition_buffer(1, :);

   % Update SOM inhibition buffer
   for buffer_i = 1 : (s2e_delay_steps-1)
       s2e_inhibition_buffer(buffer_i, :) = s2e_inhibition_buffer(buffer_i+1, :);
   end

   for buffer_i = 1 : (s2i_delay_steps-1)
       s2i_inhibition_buffer(buffer_i, :) = s2i_inhibition_buffer(buffer_i+1, :);
   end

   s2e_inhibition_buffer(s2e_delay_steps, :) = 0;
   s2i_inhibition_buffer(s2i_delay_steps, :) = 0;

   % Process SOM spikes
   if spikecount_s > 0
        if tau_rs < 0
           vs(sind_s) = vs(sind_s) - M;
        else
           vs(sind_s) = 0;
           refc_s(sind_s) = refc_s(sind_s) + tau_rs;
        end

       % Store SOM inhibition to buffer for delayed effect
       s2e_inhibition_buffer(s2e_delay_steps, :) = sum(connection_matrix_s(sind_s, 1:ne), 1);
       s2i_inhibition_buffer(s2i_delay_steps, :) = sum(connection_matrix_s(sind_s, ne+1:ne+ni), 1);
       
       % Update post-synaptic inputs for immediate SOM effects (VIP)
       hs(3,:) = hs(3,:) + sum(connection_matrix_s(sind_s, ne+ni+1:ne+ni+ns), 1);
       hv(3,:) = hv(3,:) + sum(connection_matrix_s(sind_s, ne+ni+ns+1:ne+ni+ns+nv), 1);

       % Record spike times
       spike_s(1,sind_s) = spike_s(1,sind_s) + 1;      
       spikeind = nsarray(sind_s)*spike_row + spike_s(1,sind_s) + 1;
       spike_s(round(spikeind)) = step*dt;
   end
   
   % Process VIP spikes
   if spikecount_v > 0
        if tau_rv < 0
           vv(sind_v) = vv(sind_v) - M;
        else
           vv(sind_v) = 0;
           refc_v(sind_v) = refc_v(sind_v) + tau_rv;
        end
       
       % Update post-synaptic inputs for VIP spikes
       he(3,:) = he(3,:) + sum(connection_matrix_v(sind_v, 1:ne), 1);          % VIP to E (typically minimal)
       hi(3,:) = hi(3,:) + sum(connection_matrix_v(sind_v, ne+1:ne+ni), 1);    % VIP to PV (typically minimal)
       hs(3,:) = hs(3,:) + sum(connection_matrix_v(sind_v, ne+ni+1:ne+ni+ns), 1); % VIP to SOM (important)
       hv(3,:) = hv(3,:) + sum(connection_matrix_v(sind_v, ne+ni+ns+1:ne+ni+ns+nv), 1); % VIP to VIP (minimal)

       % Record spike times
       spike_v(1,sind_v) = spike_v(1,sind_v) + 1;
       spikeind = nvarray(sind_v)*spike_row + spike_v(1,sind_v) + 1;
       spike_v(round(spikeind)) = step*dt;
   end

   % Eliminate overtime wave
   if wave_record(1) > 0 && (time - wave_record(2) > delay)   
       num_pop = sum(time - wave_record(2:wave_record(1)+1) > delay);   
       wave_record(2:wave_record(1)+1-num_pop) = wave_record(num_pop+2:wave_record(1)+1);  % eliminate overtime wave
       wave_record(1) = wave_record(1) - num_pop; 
   end

   % Start new wave
   if flag == 0 && wave_record(1) > start_threshold && time - flag_time >= flag_off_time
       flag = 1;
       flag_time = time;
       res.MFE_time(wave_count,1) = time;
       res.HEE_stat(wave_count,1) = sum(he(1,:));
       HEE_max = sum(he(1,:));
   end

   % Track maximum excitatory input
   if flag == 1 && HEE_max < sum(he(1,:))
       HEE_max = sum(he(1,:));
   end

   % End wave
   if flag == 1 && (wave_record(1) <= end_threshold || time > flag_time+10)
       flag = 0;
       flag_time = time;
       res.MFE_time(wave_count, 2) = time;
       res.HEE_stat(wave_count, 2) = HEE_max;
       res.HEE_stat(wave_count, 3) = sum(he(1,:));
       HEE_max = 0;
       wave_count = wave_count + 1;
   end
   
   % Store states in output structure
   res.VE(step,:) = ve(:);
   res.VI(step,:) = vi(:);
   res.VS(step,:) = vs(:);
   res.VV(step,:) = vv(:);
   res.HE(step,:) = [he(1,:), hi(1,:), hs(1,:), hv(1,:)];
   res.HI(step,:) = [he(2,:), hi(2,:), hs(2,:), hv(2,:)];
   res.HS(step,:) = [he(3,:), hi(3,:), hs(3,:), hv(3,:)];
   res.AS(step,:) = as(:);
   res.AV(step,:) = av(:);
   res.spikecount_e(step) = spikecount_e;
   res.spikecount_i(step) = spikecount_i;
   res.spikecount_s(step) = spikecount_s;
   res.spikecount_v(step) = spikecount_v;
end

% Assemble spike data
res.spike = [spike_e, spike_i, spike_s, spike_v];
res.spike(2:end,:) = res.spike(2:end,:)/1000;    % Convert to seconds

% Process wave data
index2 = 1;
index = 1;
while index < wave_count
   MFE_time2(index2,1) = res.MFE_time(index,1);
   HEE_stat2(index2,1) = res.HEE_stat(index,1);
   local_sp_count = 0;
   local_HEE_max = res.HEE_stat(index,2);
   while true
       local_sp_count = local_sp_count + wave_spike_count(index);
       if local_HEE_max < res.HEE_stat(index,2)
           local_HEE_max = res.HEE_stat(index,2);
       end
       HEE_stat2(index2, 3) = res.HEE_stat(index,3);
       MFE_time2(index2, 2) = res.MFE_time(index,2);
       index = index+1;
       if (res.MFE_time(index,1) - res.MFE_time(index-1,2) >= MFE_interval) || index >= wave_count
           break;
       end
   end
   wave_spike_count2(index2) = local_sp_count;
   HEE_stat2(index2,2) = local_HEE_max;
   index2 = index2 + 1;
end

wave_count = index2;
res.MFE_time = MFE_time2;
res.HEE_stat = HEE_stat2;
wave_spike_count = wave_spike_count2;
res.wave_spike_count = wave_spike_count;
res.wave_count = wave_count;

end