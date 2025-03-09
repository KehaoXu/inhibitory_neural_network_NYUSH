function real_time_visualization()
    % **Create a window (fig) - for displaying neural network simulation results and parameter adjustment**
    fig = figure('Name', 'Neural Network Simulation and Parameter Adjustment', 'Position', [50, 50, 1200, 550]);

    % Create subplot area 1 - for displaying neural network simulation results
    ax1 = subplot('Position', [0.07, 0.12, 0.5, 0.8]);
    title('Neural Network Raster Plot');

    % Store subplot handle
    handles.ax1 = ax1;
    guidata(fig, handles);  % Store GUI data

    % Create subplot area 2 - for parameter adjustment
    ax2 = subplot('Position', [0.6, 0.1, 0.4, 0.8]);
    set(ax2, 'Visible', 'off');  % Hide the axes of the subplot area

%%
    % Initialize parameters
    data = load('param2.mat'); % 'param2.mat' contains the variable 'param'
    param_som = data.param;
    param_som.duration =  1000;  % ms
    
    % Dendritic delay
    param_som.s2e_delay = 5;   % ms, from SOM to E-cells
    param_som.s2i_delay = 5;   % ms, from SOM to PV-cells
    
    % Number of Neurons
    param_som.ne = 300;
    param_som.ni = 50;
    param_som.ns = 50;
    param_som.nv = 50;     % 50 VIP neurons (similar count to SOM)
    
    % Threshold and Resting Potential
    param_som.M        = 100; % threshold potential (= 1)
    param_som.Mr       = 66;  % rest potential (= -2/3)
    
    %% refractory period
    tau = 3;
    param_som.tau_re = tau; 
    param_som.tau_ri = 2;
    param_som.tau_rs = tau;
    param_som.tau_rv = tau;  % Same refractory period as other inhibitory neurons
    
    %% probability of spike projections
    
    % CHAOS data
    % p = 0.8;        
    % param_som.p_ee = p;
    % param_som.p_ei = p;
    % param_som.p_ii = p;
    % param_som.p_ie = p;
    % param_som.p_se = p;
    % param_som.p_es = p;
    % param_som.p_is = p;
    % param_som.p_si = p;
    
    % Original Allen institute data
    % param_som.p_ee = 0.160;
    % param_som.p_ei = 0.411; 
    % param_som.p_ii = 0.451; 
    % param_som.p_ie = 0.395;
    % param_som.p_se = 0.182;
    % param_som.p_es = 0.424; 
    % param_som.p_is = 0.857; 
    % param_som.p_si = 0.030; 
    
    % Adjusted Allen institute data
    param_som.p_ee = 0.533; % 0.160
    param_som.p_ei = 0.411; % 0.411
    param_som.p_ii = 0.451; % 0.451
    param_som.p_ie = 0.220; % 0.395
    param_som.p_se = 0.182; % 0.182
    param_som.p_es = 0.424; % 0.424
    param_som.p_is = 0.857; % 0.857
    param_som.p_si = 0.030; % 0.030
    
    % Connection probabilities involving VIP
    param_som.p_sv = 0.5;   % From SOM to VIP (strong mutual inhibition)
    param_som.p_vs = 0.45;  % From VIP to SOM (strong mutual inhibition)
    param_som.p_ve = 0;     % From VIP to excitatory (minimal)
    param_som.p_vi = 0;     % From VIP to PV (minimal)
    param_som.p_ev = 0.1;   % From excitatory to VIP
    param_som.p_iv = 0;     % From PV to VIP (minimal)
    param_som.p_vv = 0;     % From VIP to VIP (minimal/zero recurrence)
    
    
    
    %% synapic strength
    external = 1; % CHAOS
    param_som.s_exe		= external;
    param_som.s_exi		= external;
    param_som.s_exs		= external;
    param_som.s_exv		= external;  % External input strength
    
    % param_som.s_ee     = 0.0200 *100  / param_som.p_ee;% 0.0180
    % param_som.s_ie     = 0.0430 *100  / param_som.p_ie;% 0.0750
    % param_som.s_ei     = 0.0240 *100  / param_som.p_ei;% 0.0240
    % param_som.s_ii     = 0.0340 *100  / param_som.p_ii;% 0.0340
    % param_som.s_si     = 0.0210 *000  / param_som.p_si;% 0.0210
    % param_som.s_es     = 0.0255 *100  / param_som.p_es;% 0.0155
    % param_som.s_se     = 0.0430 *100  / param_som.p_se;% 0.0430
    % param_som.s_is     = 0.0250 *100  / param_som.p_is;% 0.0250
     
    big = 1;
    param_som.s_ee     = 0.0180 *100 *big; % 0.0180
    param_som.s_ie     = 0.0750 *100 *big; % 0.0750
    param_som.s_ei     = 0.0240 *100 *big; % 0.0240
    param_som.s_ii     = 0.0340 *100 *big; % 0.0340
    param_som.s_si     = 0.0210 *100 *big; % 0.0210
    param_som.s_es     = 0.0155 *100 *big; % 0.0155 % 0.0455
    param_som.s_se     = 0.0430 *100 *big; % 0.0430
    param_som.s_is     = 0.0250 *100 *big; % 0.0250
    
    % param_som.s_ee     = 5.00*0.150 / param_som.p_ee;
    % param_som.s_ie     = 2.00*0.500 / param_som.p_ie;
    % param_som.s_ei     = 4.91*0.480 / param_som.p_ei;
    % param_som.s_ii     = 4.91*0.400 / param_som.p_ii;
    % % param_som.s_si     = 4.91*0.400 / param_som.p_si;
    % param_som.s_si     = 0;
    % param_som.s_es     = 4.91*0.400 / param_som.p_es;
    % param_som.s_se     = 2.00*0.500 / param_som.p_se;
    % param_som.s_is     = 4.91*0.400 / param_som.p_is;
    
    % Synaptic strengths involving VIP
    param_som.s_sv = 0.0250 *100 *big;  % From SOM to VIP (strong)
    param_som.s_vs = 0.0250 *100 *big;  % From VIP to SOM (strong)
    param_som.s_ve = 0;                 % From VIP to excitatory
    param_som.s_vi = 0;                 % From VIP to PV
    param_som.s_ev = 0.0430 *100 *big;  % From excitatory to VIP
    param_som.s_iv = 0;                 % From PV to VIP
    param_som.s_vv = 0;                 % From VIP to VIP
    
    
    %% synaptic timescale 
    % param_som.tau_ie = 2.8;   % AMPA
    % param_som.tau_ee = 2.5;    
    % param_som.tau_ei = 8.5;   % GABA
    % param_som.tau_ii = 5.8;
    % param_som.tau_se = 2.6; 
    % param_som.tau_is = 5.8;
    % param_som.tau_si = 5.8;
    % param_som.tau_es = 8.5;
    
    param_som.tau_ie = 1.2;   % AMPA
    param_som.tau_ee = 1.4;    
    param_som.tau_ei = 4.5;   % GABA
    param_som.tau_ii = 4.5;
    param_som.tau_se = 1.2; 
    param_som.tau_is = 4.5;
    param_som.tau_si = 4.5;
    param_som.tau_es = 4.5;
    
    % Time constants for VIP synapses
    param_som.tau_ev = 1.2;  % AMPA time constant
    param_som.tau_iv = 4.5;  % GABA time constant
    param_som.tau_sv = 4.5;  % GABA time constant
    param_som.tau_vv = 4.5;  % GABA time constant
    param_som.tau_ve = 1.2;  % AMPA time constant
    param_som.tau_vs = 4.5;  % GABA time constant
    param_som.tau_vi = 4.5;  % GABA time constant
    
    
    %% frequency of external input
    freq = 7000;
    param_som.lambda_e = freq;    
    param_som.lambda_i = freq;
    param_som.lambda_s = freq;
    param_som.lambda_v = freq;  % External input frequency
    
    
    %% SFA parameters for VIP and SOM neurons
    param_som.tau_av = 50;  % Adaptation time constant (ms)
    param_som.b_v = 0.8;    % Adaptation strength 
    param_som.tau_as = 50;  % Adaptation time constant (ms)
    param_som.b_s = 0.8;    % Adaptation strength 

    
    %%
    % **Create slider controls**
    label_pos = [0.60, 0.87, 0.07, 0.05];
    edit_pos = [0.67, 0.87, 0.05, 0.05];
    create_slider(fig, 'duration', label_pos, edit_pos, 500, 5000, param_som.duration, @update_duration);
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 's2e_delay', label_pos, edit_pos, 1, 10, param_som.s2e_delay, @update_s2e_delay);
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 's2i_delay', label_pos, edit_pos, 1, 10, param_som.s2i_delay, @update_s2i_delay);
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'ne', label_pos, edit_pos, 100, 1000, param_som.ne, @update_ne);
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'ni', label_pos, edit_pos, 10, 500, param_som.ni, @update_ni);
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'ns', label_pos, edit_pos, 10, 500, param_som.ns, @update_ns);
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'nv', label_pos, edit_pos, 10, 500, param_som.nv, @update_nv);  % New parameter nv
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'M', label_pos, edit_pos, 0, 100, param_som.M, @update_M);
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'Mr', label_pos, edit_pos, 0, 100, param_som.Mr, @update_Mr);

    label_pos = [0.82, 0.87, 0.05, 0.05];
    edit_pos = [0.87, 0.87, 0.05, 0.05];
    create_slider(fig, 's_es', label_pos, edit_pos, 0, 100, param_som.s_es, @update_s_es);  % New parameter s_es
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 's_is', label_pos, edit_pos, 0, 100, param_som.s_is, @update_s_is);  % New parameter s_is
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 's_sv', label_pos, edit_pos, 0, 100, param_som.s_sv, @update_s_sv);  % New parameter s_sv
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 's_vs', label_pos, edit_pos, 0, 100, param_som.s_vs, @update_s_vs);  % New parameter s_vs
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'tau_av', label_pos, edit_pos, 0, 100, param_som.tau_av, @update_tau_av);  % New parameter tau_av
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'tau_as', label_pos, edit_pos, 0, 100, param_som.tau_as, @update_tau_as);  % New parameter tau_as
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'b_v', label_pos, edit_pos, 0, 1, param_som.b_v, @update_b_v);  % New parameter b_v
    label_pos = label_pos - [0, 0.1, 0, 0];
    edit_pos = edit_pos - [0, 0.1, 0, 0];
    create_slider(fig, 'b_s', label_pos, edit_pos, 0, 1, param_som.b_s, @update_b_s);  % New parameter b_s

    % **Run simulation button**
    uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Run Simulation', ...
              'Units', 'normalized', 'Position', [0.85, 0.07, 0.12, 0.05], 'FontSize', 12, ...
              'Callback', @run_simulation);

    % **Run simulation**
    function run_simulation(~, ~)
        handles = guidata(fig);
        cla(handles.ax1);

        tic;
        res_lif = model_LIF_SOM_VIP(param_som,[]);
        toc;

        % **Switch to fig to plot the image**
        figure(fig);
        
        subplot(handles.ax1);
        rasterplot3(res_lif, param_som);
        title('Neural Network Raster Plot');
    end

    % **Parameter update functions**
    function update_duration(hObj, ~)
        param_som.duration = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_s2e_delay(hObj, ~)
        param_som.s2e_delay = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_s2i_delay(hObj, ~)
        param_som.s2i_delay = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_ne(hObj, ~)
        param_som.ne = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_ni(hObj, ~)
        param_som.ni = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_ns(hObj, ~)
        param_som.ns = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_nv(hObj, ~)
        param_som.nv = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_M(hObj, ~)
        param_som.M = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_Mr(hObj, ~)
        param_som.Mr = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_s_es(hObj, ~)
        param_som.s_es = get(hObj, 'Value');
        update_slider_label(hObj);
    end

    function update_s_is(hObj, ~)
        param_som.s_is = get(hObj, 'Value');
        update_slider_label(hObj);
    end

    function update_s_sv(hObj, ~)
        param_som.s_sv = get(hObj, 'Value');
        update_slider_label(hObj);
    end

    function update_s_vs(hObj, ~)
        param_som.s_vs = get(hObj, 'Value');
        update_slider_label(hObj);
    end

    function update_tau_av(hObj, ~)
        param_som.tau_av = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_tau_as(hObj, ~)
        param_som.tau_as = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_b_v(hObj, ~)
        param_som.b_v = get(hObj, 'Value');
        update_slider_label(hObj);
    end

    function update_b_s(hObj, ~)
        param_som.b_s = get(hObj, 'Value');
        update_slider_label(hObj);
    end

    function update_slider_label(hObj)
        label = get(hObj, 'UserData');
        set(label, 'String', sprintf('%s', get(label, 'UserData')));
    end
end

% **Modify `create_slider` to place the slider to the right of the text and add an input box**
function create_slider(parent, label, labelPos, editBoxPos, minVal, maxVal, initVal, callback)
    % **Label (left-aligned)**
    labelHandle = uicontrol('Parent', parent, 'Style', 'text', 'Units', 'normalized', 'Position', labelPos, ...
              'String', sprintf('%s', label), 'FontSize', 12, 'HorizontalAlignment', 'left', ...
              'UserData', label);

    % **Slider (right side)**
    sliderPos = [editBoxPos(1) + editBoxPos(3) + 0.01, editBoxPos(2), 0.05, editBoxPos(4)];
    slider = uicontrol('Parent', parent, 'Style', 'slider', 'Units', 'normalized', 'Min', minVal, 'Max', maxVal, ...
              'Value', initVal, 'Position', sliderPos, 'Callback', callback, 'UserData', labelHandle);

    % **Input box (middle)**
    editBox = uicontrol('Parent', parent, 'Style', 'edit', 'Units', 'normalized', 'Position', editBoxPos, ...
              'String', num2str(initVal), 'FontSize', 10, 'Callback', @(hObj, event) edit_callback(hObj, slider));

    % **Listen to slider value changes and update the label and input box**
    addlistener(slider, 'ContinuousValueChange', @(hObj, event) ...
        set(labelHandle, 'String', sprintf('%s', label)));
    addlistener(slider, 'ContinuousValueChange', @(hObj, event) ...
        set(editBox, 'String', num2str(round(get(hObj, 'Value')))));

    % **Input box callback function**
    function edit_callback(hObj, slider)
        val = str2double(get(hObj, 'String'));
        if isnan(val) || val < minVal || val > maxVal
            set(hObj, 'String', num2str(get(slider, 'Value')));
        else
            set(slider, 'Value', val);
            callback(slider, []);
        end
    end
end
