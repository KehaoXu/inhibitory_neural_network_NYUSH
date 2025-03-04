function real_time_visualization()
    % **创建一个窗口 (fig) - 用于显示神经网络仿真结果和参数调节**
    fig = figure('Name', '神经网络仿真与参数调节', 'Position', [50, 50, 1200, 550]);

    % 创建子图区域1 - 用于显示神经网络仿真结果
    ax1 = subplot('Position', [0.07, 0.12, 0.6, 0.8]);
    title('神经网络光栅图');

    % 存储 subplot 句柄
    handles.ax1 = ax1;
    guidata(fig, handles);  % 存储 GUI 数据

    % 创建子图区域2 - 用于参数调节
    ax2 = subplot('Position', [0.7, 0.1, 0.3, 0.8]);
    set(ax2, 'Visible', 'off');  % 隐藏子图区域的坐标轴

%%
    % 初始化参数
    data = load('param2.mat'); % 'param2.mat' 中包含变量 'param'
    param_som = data.param;
    param_som.duration =  1000;  % ms
    
    param_som.s2e_delay = 5;   % ms, from SOM to E-cells
    param_som.s2i_delay = 5;   % ms, from SOM to PV-cells
    
    param_som.ne = 300;
    param_som.ni = 50;
    param_som.ns = 50;
    
    param_som.M        = 100; % threshold potential (= 1)
    param_som.Mr       = 66;  % rest potential (= -2/3)
    
    %% refractory period
    tau = 3;
    param_som.tau_re = tau; 
    param_som.tau_ri = 2;
    param_som.tau_rs = tau;
    
    %% probability of spike projections
    
    % Original Allen institute data
    param_som.p_ee = 0.160;
    param_som.p_ei = 0.411; 
    param_som.p_ii = 0.451; 
    param_som.p_ie = 0.395;
    param_som.p_se = 0.182;
    param_som.p_es = 0.424; 
    param_som.p_is = 0.857; 
    param_som.p_si = 0.030; 
    
    %% synapic strength
    external = 1; % CHAOS
    param_som.s_exe    = external;
    param_som.s_exi    = external;
    param_som.s_exs    = external;
     
    big = 1;
    param_som.s_ee     = 0.0180 *100 *big; % 0.0180
    param_som.s_ie     = 0.0750 *100 *big; % 0.0750
    param_som.s_ei     = 0.0240 *100 *big; % 0.0240
    param_som.s_ii     = 0.0340 *100 *big; % 0.0340
    param_som.s_si     = 0.0210 *100 *big; % 0.0210
    param_som.s_es     = 0.0155 *100 *big; % 0.0155 % 0.0455
    param_som.s_se     = 0.0430 *100 *big; % 0.0430
    param_som.s_is     = 0.0250 *100 *big; % 0.0250
    
    %% synaptic timescale 
    param_som.tau_ie = 1.2;   % AMPA
    param_som.tau_ee = 1.4;    
    param_som.tau_ei = 4.5;   % GABA
    param_som.tau_ii = 4.5;
    param_som.tau_se = 1.2; 
    param_som.tau_is = 4.5;
    param_som.tau_si = 4.5;
    param_som.tau_es = 4.5;
    
    %% frequency of exteranl input
    freq = 7000;
    param_som.lambda_e = freq;    
    param_som.lambda_i = freq;
    param_som.lambda_s = freq;

    % **创建滑块控件**
    create_slider(fig, 'duration', [0.70, 0.85, 0.1, 0.05], [0.80, 0.85, 0.1, 0.05], 500, 5000, param_som.duration, @update_duration);
    create_slider(fig, 's2e_delay', [0.70, 0.75, 0.1, 0.05], [0.80, 0.75, 0.1, 0.05], 1, 10, param_som.s2e_delay, @update_s2e_delay);
    create_slider(fig, 's2i_delay', [0.70, 0.65, 0.1, 0.05], [0.80, 0.65, 0.1, 0.05], 1, 10, param_som.s2i_delay, @update_s2i_delay);
    create_slider(fig, 'ne', [0.70, 0.55, 0.1, 0.05], [0.80, 0.55, 0.1, 0.05], 100, 1000, param_som.ne, @update_ne);
    create_slider(fig, 'ni', [0.70, 0.45, 0.1, 0.05], [0.80, 0.45, 0.1, 0.05], 10, 500, param_som.ni, @update_ni);
    create_slider(fig, 'ns', [0.70, 0.35, 0.1, 0.05], [0.80, 0.35, 0.1, 0.05], 10, 500, param_som.ns, @update_ns);
    create_slider(fig, 'M', [0.70, 0.25, 0.1, 0.05], [0.80, 0.25, 0.1, 0.05], 0, 100, param_som.M, @update_M);
    create_slider(fig, 'Mr', [0.70, 0.15, 0.1, 0.05], [0.80, 0.15, 0.1, 0.05], 0, 100, param_som.Mr, @update_Mr);

    % **运行仿真的按钮**
    uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', '运行仿真', ...
              'Units', 'normalized', 'Position', [0.80, 0.05, 0.1, 0.05], 'FontSize', 12, ...
              'Callback', @run_simulation);

    % **运行仿真**
    function run_simulation(~, ~)
        handles = guidata(fig);
        cla(handles.ax1);

        tic;
        res_lif = model_LIF_SOM(param_som, []);
        toc;

        % **切换到 fig 绘制图像**
        figure(fig);
        
        subplot(handles.ax1);
        rasterplot2(res_lif, param_som);
        title('神经网络光栅图');
    end

    % **参数更新函数**
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

    function update_M(hObj, ~)
        param_som.M = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_Mr(hObj, ~)
        param_som.Mr = round(get(hObj, 'Value'));
        update_slider_label(hObj);
    end

    function update_slider_label(hObj)
        label = get(hObj, 'UserData');
        set(label, 'String', sprintf('%s: %d', get(label, 'UserData'), round(get(hObj, 'Value'))));
    end
end

% **修改 `create_slider` 使滑块在文字右侧，并添加输入框**s
function create_slider(parent, label, labelPos, sliderPos, minVal, maxVal, initVal, callback)
    % **标签（靠左）**
    labelHandle = uicontrol('Parent', parent, 'Style', 'text', 'Units', 'normalized', 'Position', labelPos, ...
              'String', sprintf('%s: %d', label, initVal), 'FontSize', 12, 'HorizontalAlignment', 'left', ...
              'UserData', label);

    % **滑块（中间）**
    slider = uicontrol('Parent', parent, 'Style', 'slider', 'Units', 'normalized', 'Min', minVal, 'Max', maxVal, ...
              'Value', initVal, 'Position', sliderPos, 'Callback', callback, 'UserData', labelHandle);
    
    % % **调整滑块拇指大小**
    % jSlider = findjobj(slider);
    % jSlider.setThumbSize(20);  % 设置滑块拇指大小

    % **输入框（右侧）**
    editBoxPos = [sliderPos(1) + sliderPos(3) + 0.02, sliderPos(2), 0.05, sliderPos(4)];
    editBox = uicontrol('Parent', parent, 'Style', 'edit', 'Units', 'normalized', 'Position', editBoxPos, ...
              'String', num2str(initVal), 'FontSize', 10, 'Callback', @(hObj, event) edit_callback(hObj, slider));

    % **监听滑块值的变化，更新标签和输入框**
    addlistener(slider, 'ContinuousValueChange', @(hObj, event) ...
        set(labelHandle, 'String', sprintf('%s: %d', label, round(get(hObj, 'Value')))));
    addlistener(slider, 'ContinuousValueChange', @(hObj, event) ...
        set(editBox, 'String', num2str(round(get(hObj, 'Value')))));

    % **输入框回调函数**
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
