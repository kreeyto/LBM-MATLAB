% Roofline Model - NVIDIA RTX 4090 (MATLAB version)

% GPU specs
peak_flops = 82.58;         % TFLOP/s
peak_bandwidth = 1.008;     % TB/s

% Operational intensity range (FLOP/byte)
oi = logspace(-1, 4, 500);
roofline = min(peak_flops, peak_bandwidth .* oi);

% User's data point
user_oi = 2206 / 520;       % FLOP/byte
user_flops = 3.31;          % TFLOP/s

% Montessori points
labels = {'CUDA 2 fluids', 'OACC 2 fluids', 'CUDA single fluid', 'OACC single fluid'};
oi_mont = [0.9, 1.7, 2.0, 2.2];
flops_mont = [0.65, 1.3, 1.5, 1.7];
colors = {'red', 'green', 'yellow', 'magenta'};

% Plot
figure;
loglog(oi, roofline, 'r-', 'LineWidth', 2); hold on;
loglog(user_oi, user_flops, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
text(user_oi*1.1, user_flops*1.1, 'This Work', 'Color', 'b');

% Montessori points
for i = 1:length(labels)
    loglog(oi_mont(i), flops_mont(i), 'o', 'MarkerFaceColor', colors{i});
    text(oi_mont(i)*1.1, flops_mont(i)*1.1, labels{i}, 'FontSize', 9);
end

% Labels and title
xlabel('Operational Intensity (FLOP/byte)');
ylabel('Performance (TFLOP/s)');
title('Roofline Model - NVIDIA RTX 4090');
grid on;
legend('RTX 4090 Roofline', 'This Work', 'Location', 'northwest');
xlim([0.1, 1e4]);
ylim([0.1, 100]);
set(gca, 'XScale', 'log', 'YScale', 'log');
