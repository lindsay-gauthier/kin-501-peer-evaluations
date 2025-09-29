%Create signals using the mathematical formula x(t) = A • cos(2πf0t + φ) 

t = 0:0.001:6; % 6s time with time steps of 0.001s 
% Manipulate these variables to generate new signals using 'for loops'


% 125 graphs in one figure?
% plot the effect in ONE graph?
% define the signal:EMG in mV; Arbitrary unit
% axis increment
% for A,f0,phase; [inside loops] save index array; end end end; plot outside
% have array & plot later
% animation? show one on top of others
% legend outside
% instead of color, use shading to show the effect?
% change unit for the graph & change title??
% 3min presentation


%% create all 125 plots
TIME = 0:0.001:6
AMPLITUDE_VALUES = [1, 2, 3, 4, 5];
FREQUENCY_VALUES = [0.2, 0.5, 1, 2, 3];
% 0.1-0.5Hz for natural CoP sway in quiet standing (ankle strategies)
% 1-3Hz for CoP sway in more chanllenging standing (hip strategies)
% 5Hz for EMG activity - reflexive neuromuscular control (stretch reflex)
PHASE_VALUES = [0, 1/4*pi, 1/2*pi, 3/4*pi, pi];
results = cell(length(AMPLITUDE_VALUES)*length(FREQUENCY_VALUES)*length(PHASE_VALUES), 4);
row = 1

for iAMPLITUDE = AMPLITUDE_VALUES;
    for iFREQUENCY = FREQUENCY_VALUES;
        for iPHASE = PHASE_VALUES;
            x = iAMPLITUDE * cos(2*pi*iFREQUENCY*TIME+ iPHASE);
            results{row,1} = x;              
            results{row,2} = iAMPLITUDE;
            results{row,3} = iFREQUENCY;
            results{row,4} = iPHASE;
            row = row + 1;
        end
    end
end

%% 3D plots of variable effect
% --- Parameters ---
% X axis: TIME
% Y levels: AMPLITUDE_VALUES
% FREQUENCY_VALUES
phi1 = 0;                                  % fix phase 
phi2 = pi;
% --- Build surface grids ---
[T, Agrid] = meshgrid(TIME, AMPLITUDE_VALUES);           % T: [numA x numT], Agrid same size

% Row-wise frequency map (same length as A)
FREQUENCY_LIST = FREQUENCY_VALUES(:);           % column vector [numA x 1]      
Fgrid = repmat(FREQUENCY_LIST, 1, 6001)         % expand to [numA x numT]
                                                % 6001 = numel(TIME)

%PHASE_LIST = PHASE_VALUES(:);               
%Pgrid = repmat(PHASE_LIST, 1, 6001)
% --- Signal surface (Z) ---
SIGNAL1 = Agrid .* cos(2*pi.*Fgrid.*T + phi1);   


% --- Greyscale coloring: dark at peaks, pale near zero ---
C = abs(SIGNAL1);                          % intensity by magnitude
C = C ./ max(C(:));                        % normalize 0.1
C = 1 - C;                                 % invert so peaks -> dark, zero -> pale

% --- Plot as smooth grey surf ---
figure('Color','w');

surf(T, Agrid, SIGNAL1,C, 'EdgeColor','none');  % use C as color data (greyscale)
colormap(gray);                            % single grey palette
shading interp;                            % smooth color transitions

                   

view([-35 30]);                            % nice 3D angle suggested
box on; grid on;


xlabel('Time (s)');
ylabel('Amplitude (a.u.)');                
zlabel('Signal');
title('The Effect of Amplitude and Frequency on CoP Sway Cosine Waves', fontsize = 20);

% Optional subtle lighting (keeps greyscale look)
% camlight headlight; material dull;

%%
% --- Parameters ---
phi1 = 0; 
phi2 = pi;

% --- Build surface grids ---
[T, Agrid] = meshgrid(TIME, AMPLITUDE_VALUES);            % [nA x nT]
FREQUENCY_LIST = FREQUENCY_VALUES(:);                     % [nA x 1]
Fgrid = repmat(FREQUENCY_LIST, 1, numel(TIME));           % [nA x nT]

% --- Signals ---
SIGNAL1 = Agrid .* cos(2*pi.*Fgrid.*T + phi1);
SIGNAL2 = Agrid .* cos(2*pi.*Fgrid.*T + phi2);

% --- Greyscale coloring for base surface (pale near zero, dark at peaks) ---
C1 = abs(SIGNAL1);
C1 = C1 ./ max(C1(:));
C1 = 1 - C1;

% --- Plot ---
F = figure('Color','w'); hold on
s1 = surf(T, Agrid, SIGNAL1, C1, 'EdgeColor','none');     % base (greyscale)
colormap(gray);                                           % single colormap
shading interp

% --- Red overlay: solid red, transparency ~ |signal| ---
A2 = abs(SIGNAL2);                    
A2 = A2 ./ max(A2(:));                % 0..1
gamma = 1.0;                          % try 1.5 to make mid-values paler
A2 = A2.^gamma;

s2 = surf(T, Agrid, SIGNAL2, 'EdgeColor','none', ...
          'FaceColor',[1 0 0], ...    % fixed red (doesn't use colormap)
          'FaceAlpha','interp', ...
          'AlphaData', A2);            % opaque at peaks, transparent near zero

view([-35 30]); box on; grid on
xlabel('Time (s)', FontSize=15)
ylabel('Amplitude (a.u.)', FontSize=15)    
zlabel('Signal')
title('The Effect of Amplitude, Frequency, and Phase on COP Sway Cosine Waves', fontsize = 20)
camlight headlight; material dull
legend('Phase = 0', 'Phase = π', 'box', 'off', FontSize=10)

%% 2D plots of variable effect
figure;

subplot(3, 1, 1);
plot(TIME, results{1,1}, LineWidth=1.5);
hold on;
plot(TIME, results{2,1}, LineWidth=1.5);
plot(TIME, results{3,1}, LineWidth=1.5);
plot(TIME, results{4,1}, LineWidth=1.5);
plot(TIME, results{5,1}, LineWidth=1.5);
box off;
title('Effect of Phase on Signal', fontsize = 24);
legend("phase = 0", "phase = 1/2π", "phase = 1/4π", "phase = 3/4π","phase = π", 'Box','off', 'Location','eastoutside');
xlabel('Time (s)', fontsize = 15);
ylabel('CoP Displacement_A_P (cm)', fontsize = 15);
subplot(3, 1, 2);
plot(TIME, results{1,1}, LineWidth=1.5);
hold on;
plot(TIME, results{6,1}, LineWidth=1.5);
plot(TIME, results{11,1}, LineWidth=1.5);
plot(TIME, results{16,1}, LineWidth=1.5);
plot(TIME, results{21,1}, LineWidth=1.5);
box off;
title('Effect of Frequency (f0) on Signal', fontsize = 24);
legend("f0 = 0.2Hz", "f0 = 0.5Hz", "f0 = 1Hz", "f0 = 2Hz","f0 = 3Hz", 'Box','off', 'Location','eastoutside');
xlabel('Time (s)', fontsize = 15);
ylabel('CoP Displacement_A_P (cm)', fontsize = 15);
subplot(3, 1, 3);
plot(TIME, results{1,1}, LineWidth=1.5);
hold on;
plot(TIME, results{26,1}, LineWidth=1.5);
plot(TIME, results{51,1}, LineWidth=1.5);
plot(TIME, results{76,1}, LineWidth=1.5);
plot(TIME, results{101,1}, LineWidth=1.5);
box off;
title('Effect of Amplitude (A) on Signal', fontsize = 24);
legend("A = 1 a.u.", "A = 2", "A = 3", "A = 4","A = 5", 'Box','off', 'Location','eastoutside');
xlabel('Time (s)', fontsize = 15);
ylabel('CoP Displacement_A_P (cm)', fontsize = 15);
%%
tol = 1e-6;
order1 = find(abs(results{1,1}) < tol);
zerovalue1 = TIME(order1(1,1));
plot(zerovalue1, results{1,1}(order1), 'ro');
order2 = find(abs(results{2,1}) < tol);
zerovalue2 = TIME(order2(1,2));
plot(zerovalue2, results{2,1}(order2), 'ro');
order3 = find(abs(results{3,1}) < tol);
zerovalue3 = TIME(order3(1,3));
plot(zerovalue3, results{3,1}(order3), 'ro');
order4 = find(abs(results{4,1}) < tol);
zerovalue4 = TIME(order4(1,4));
plot(zerovalue4, results{4,1}(order4), 'ro');
order5 = find(abs(results{5,1}) < tol);
zerovalue5 = TIME(order5(1,5));
plot(zerovalue5, results{5,1}(order5), 'ro');
box off
%legend("phase = 0", 'Box','off', 'Location','eastoutside');
hold off;
subplot(5, 1, 2);
plot(TIME, results{2,1});
hold on;
order2 = find(abs(results{2,1}) < tol);
zerovalue2 = TIME(order2(1,2));
plot(zerovalue2, results{2,1}(order2), 'ro', 'MarkerFaceColor','r');
box off
%legend("phase = 1/2π", 'Box','off', 'Location','eastoutside');
hold off;
subplot(5, 1, 3);
plot(TIME, results{3,1});
hold on;
order3 = find(abs(results{3,1}) < tol);
zerovalue3 = TIME(order3(1,3));
plot(zerovalue3, results{3,1}(order3), 'ro', 'MarkerFaceColor','r');
box off
%legend("phase = 1/4π", 'Box','off', 'Location','eastoutside');
hold off;
subplot(5, 1, 4);
plot(TIME, results{4,1});
hold on;
order4 = find(abs(results{4,1}) < tol);
zerovalue4 = TIME(order4(1,4));
plot(zerovalue4, results{4,1}(order4), 'ro', 'MarkerFaceColor','r');
box off
%legend("phase = 3/4π", 'Box','off', 'Location','eastoutside');
hold off;
subplot(5, 1, 5);
plot(TIME, results{5,1});
hold on;
order5 = find(abs(results{5,1}) < tol);
zerovalue5 = TIME(order5(1,5));
plot(zerovalue5, results{5,1}(order5), 'ro', 'MarkerFaceColor','r');
box off
%legend("phase = π", 'Box','off', 'Location','eastoutside');
xlabel('Time (s)', fontsize = 15);
hold off;
X = axes(fig,'visible','off'); 
X.YLabel.Visible = 'on'; 
ylabel(X, 'CoP Displacement (cm)', fontsize = 15); 
%%
a = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
b = [1,1,0,0,1,1,0,0,1,1,0,0,1,1];
plot(a, b, 'LineWidth', 4)


%%
%plot(TIME, results{36,1}, 'k', LineWidth=4);
a =[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
b=[1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0]
plot(a,b, LineWidth=4)
%xlim([0,2000])
%hold on
%plot(TIME, results{2,1}, 'b', LineWidth=1.5);
%xlabel('Time (s)', fontsize = 15, FontWeight='bold');
%ylabel('CoP Displacement_A_P (cm)', fontsize = 15, FontWeight='bold');
%set(gca, 'LineWidth', 2);   % 2 is a good starting point
axis off
box off
%hold off

%% Effect of Amplitude 2D
f0 = 1;
phase = 0;
A_values = [1, 2, 3, 4, 5];
figure;
hold on;
for A = A_values;
    x = A * cos(2*pi*f0*t+ phase)
    plot (t,x, linewidth=1.5)
end
legend('A = 1', 'A = 2', 'A = 3', 'A = 4', 'A = 5');
xlabel('Time (s)');
ylabel('x(t)');
title('Effect of Amplitude (A) on Signal');
hold off
%% Effect of Frequency 2D
A = 2;
phase = 0;
f0_values = [0.1, 0.2, 0.5, 1, 2];
figure;
hold on;
for f0 = f0_values;
    x = A * cos(2*pi*f0*t+ phase)
    plot (t,x, linewidth=1.5)
end
legend('f0 = 0.1', 'f0 = 0.2', 'f0 = 0.5', 'f0 = 1', 'f0 = 2');
xlabel('Time (s)');
ylabel('x(t)');
title('Effect of frequency (f0) on Signal');
hold off
%% Effect of Phase 2D
A = 2; 
f0 = 1;
phase_values = [0, 1/2*pi, 1/4*pi, 3/4*pi, pi];
figure;
hold on;
for phase = phase_values;
    x = A * cos(2*pi*f0*t+phase)
    plot (t,x, linewidth=1.5)
end
legend('f0 = 0.1', 'f0 = 0.2', 'f0 = 0.5', 'f0 = 1', 'f0 = 2');
xlabel('Time (s)');
ylabel('x(t)');
title('Effect of phase on Signal');
hold off
%% Amplitude Effect 3D
t = 0:0.001:6;    % time
f0 = 2; phase = 0;
A1 = 1

y = A1 * ones(size(t)); % x-axis = constant amplitude
x = t;  % y-axis = time
z = A1 * cos(2*pi*f0*t + phase);   % z-axis = signal value

figure;
plot3(x, y, z, 'LineWidth', 1.5)
hold on;

A2 = 2
y = A2 * ones(size(t)); % y-axis = constant amplitude
z = A2 * cos(2*pi*f0*t + phase);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

A3 = 3
y = A3 * ones(size(t)); % y-axis = constant amplitude
z = A3 * cos(2*pi*f0*t + phase);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

A4 = 4
y = A4 * ones(size(t)); % y-axis = constant amplitude
z = A4 * cos(2*pi*f0*t + phase);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

A5 = 5
y = A5 * ones(size(t)); % y-axis = constant amplitude
z = A5 * cos(2*pi*f0*t + phase);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

grid on
xlabel('Time (s)')
ylabel('Amplitude (A)')
zlabel('Signal x(t)')
title('Effect of Amplitude (A) on Signal')
hold off

%% Frequency Effect 3D
t = 0:0.001:6;    % time
A = 2; phase = 0;
f01 = 0.1

y = f01 * ones(size(t)); % x-axis = constant frequency value
x = t;  % y-axis = time
z = A * cos(2*pi*f01*t + phase);   % z-axis = signal value

figure;
plot3(x, y, z, 'LineWidth', 1.5)
hold on;

f02 = 0.2
y = f02 * ones(size(t)); % y-axis = constant amplitude
z = A * cos(2*pi*f02*t + phase);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

f03 = 0.5
y = f03 * ones(size(t)); % y-axis = constant amplitude
z = A * cos(2*pi*f03*t + phase);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

f04 = 1
y = f04 * ones(size(t)); % y-axis = constant amplitude
z = A * cos(2*pi*f04*t + phase);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

f05 = 2
y = f05 * ones(size(t)); % y-axis = constant amplitude
z = A * cos(2*pi*f05*t + phase);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

grid on
xlabel('Time (s)')
ylabel('Frequency (f0)')
zlabel('Signal x(t)')
title('Effect of Frequency (f0) on Signal')
hold off

%% Phase Effect 3D
t = 0:0.001:6;    % time
A = 2; f0 = 0.5;

phase1 = 0
y = phase1 * ones(size(t)); % x-axis = constant frequency value
x = t;  % y-axis = time
z = A * cos(2*pi*f0*t + phase1);   % z-axis = signal value

figure;
plot3(x, y, z, 'LineWidth', 1.5)
hold on;

phase2 = 1/2*pi;
y = phase2 * ones(size(t)); % y-axis = constant amplitude
z = A * cos(2*pi*f0*t + phase2);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

phase3 = 1/4*pi;
y = phase3 * ones(size(t)); % y-axis = constant amplitude
z = A * cos(2*pi*f0*t + phase3);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

phase4 = 3/4*pi;
y = phase4 * ones(size(t)); % y-axis = constant amplitude
z = A * cos(2*pi*f0*t + phase4);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

phase5 = pi;
y = phase5 * ones(size(t)); % y-axis = constant amplitude
z = A * cos(2*pi*f0*t + phase5);   % z-axis = new signal value
plot3(x,y,z, 'LineWidth', 1.5);

grid on
xlabel('Time (s)')
ylabel('Phase')
zlabel('Signal x(t)')
title('Effect of Phase on Signal')
hold off