%% eegmovie_testing.m

%% 2-D Scalp Topography Animation using eegmovie (with frame rate control and MP4 export)

% Load your data
load('/Volumes/T7/ERP Files/Epoched Files/paired_ttest/LIMO.mat');
data = load('/Volumes/T7/ERP Files/Epoched Files/paired_ttest/paired_samples_ttest_parameter_14.mat');
paired_samples = data.paired_samples;
average_files = {'/Volumes/T7/ERP Files/Epoched Files/derivatives/parameter_1_Mean_of_Mean.mat', '/Volumes/T7/ERP Files/Epoched Files/derivatives/parameter_4_Mean_of_Mean.mat'};

% Load the two Mean_of_Mean files to get the voltage difference
data1_struct = load(average_files{1});
data2_struct = load(average_files{2});

% Extract the average voltage values (index 2 of the last dimension) and compute difference
voltage_data1 = squeeze(data1_struct.Data.mean(:,:,:,2));
voltage_data2 = squeeze(data2_struct.Data.mean(:,:,:,2));
difference_voltage = voltage_data1 - voltage_data2;

% Extract t-values (4th dimension)
plot_data = squeeze(paired_samples(:,:,4)); % t-values

% Set up time parameters
times = linspace(LIMO.data.start, LIMO.data.end, size(plot_data,2));

% Determine color limits (symmetric around zero for t-values)
abs_max = max(abs(plot_data(:)));
colorbar_limits = [-abs_max, abs_max];

% Determine color limits for the voltage difference plot
volt_abs_max = max(abs(difference_voltage(:)));
volt_colorbar_limits = [-volt_abs_max, volt_abs_max];

% Calculate the actual sampling rate from your data
actual_duration_ms = LIMO.data.end - LIMO.data.start; % Duration in ms
actual_duration_sec = actual_duration_ms / 1000; % Convert to seconds
num_samples = size(plot_data, 2);
actual_srate = num_samples / actual_duration_sec; % Samples per second

% ===== PLAYBACK PARAMETERS =====
% Time multiplier for playback speed (0.5 = half speed, 1 = real time, 2 = double speed)
time_multiplier = 0.0625; % 0.0625

% Target frame rate for the output video (e.g., 30 fps is standard for video)
% Lower values = fewer frames to generate (faster), but less smooth
% Higher values = more frames to generate (slower), but smoother
target_frame_rate = 50; % 50

% Calculate which frames to extract from the data
% We need target_frame_rate frames per second for the output duration
output_duration_sec = actual_duration_sec / time_multiplier;
total_frames_needed = round(target_frame_rate * output_duration_sec);

% Calculate which data indices to use for each frame
frame_indices = round(linspace(1, num_samples, total_frames_needed));
frame_indices = unique(frame_indices); % Remove any duplicates

% Update total frames needed after removing duplicates
total_frames_needed = length(frame_indices);

% Define which frames to animate
movieframes = frame_indices; % Selected frames based on target frame rate

% Display frame generation info
fprintf('\n=== Frame Generation Settings ===\n');
fprintf('Original data points: %d\n', num_samples);
fprintf('Target frame rate: %d fps\n', target_frame_rate);
fprintf('Frames to generate: %d\n', total_frames_needed);
fprintf('Frame reduction: %.1f%% (%.1fx faster generation)\n', ...
   (1 - total_frames_needed/num_samples) * 100, num_samples/total_frames_needed);
fprintf('Estimated generation time: ~%.1f seconds\n', total_frames_needed / 10); % Assuming ~10 fps generation
fprintf('==================================\n\n');

% Create the movie frames using eegmovie
figure('Position', [100 100 800 600]);

% LIMO-style colormap
% Colormap for the second plot (t-values)
t_value_colormap = limo_color_images(plot_data(~isnan(plot_data)));

% Call eegmovie with correct parameters
[Movie, Colormap] = eegmovie_parallel_dual(difference_voltage, actual_srate, LIMO.data.chanlocs, ...
    'title', sprintf('Paired T-test Animation (%.1fx speed, %d fps)', time_multiplier, target_frame_rate), ...
    'data2', plot_data, ...
    'mode', '2D', ...
    'layout', 'dual', ...
    'title', 'Voltage Difference & Test Statistic', ...
    'subtitle1', 'Voltage Difference (μV)', ...
    'subtitle2', 'T-values', ...
    'colormap2', t_value_colormap, ...
    'minmax', volt_colorbar_limits, ...
    'minmax2', colorbar_limits, ...
    'startsec', LIMO.data.start/1000, ...
    'movieframes', frame_indices, ...
    'headlinewidth', 12, ...
    'topo_linewidth', 3, ...
    'timecourse', 'off', ...
    'framenum', 'off', ...
    'time', 'on', ...
    'showcolorbars', 'on', ...
    'resolution', [2000 1000]);

% Apply colormap
colormap(jet);

% Display information
fprintf('\n=== Movie Information ===\n');
fprintf('Animation created with %d frames\n', size(Movie,2));
fprintf('Data duration: %.1f ms (%.2f seconds)\n', actual_duration_ms, actual_duration_sec);
fprintf('Time range: %.1f to %.1f ms\n', LIMO.data.start, LIMO.data.end);
fprintf('Playback speed multiplier: %.1fx\n', time_multiplier);
fprintf('Output frame rate: %d fps\n', target_frame_rate);
fprintf('Output duration: %.2f seconds\n', output_duration_sec);
fprintf('Color scale: %.2f to %.2f (t-values)\n', colorbar_limits(1), colorbar_limits(2));
fprintf('=========================\n\n');

% Save the movie as MP4 with proper frame rate
output_filename = sprintf('paired_ttest_animation_%.1fx_speed_%dfps.mp4', time_multiplier, target_frame_rate);

% Create VideoWriter object for MP4
videoObj = VideoWriter(output_filename, 'MPEG-4');
videoObj.FrameRate = target_frame_rate; % Set the target frame rate
videoObj.Quality = 95; % High quality (1-100 scale)

% Open the video file for writing
open(videoObj);

% Write each frame to the video
fprintf('Writing MP4 file...\n');
for i = 1:size(Movie, 2)
   writeVideo(videoObj, Movie(:,i));
   if mod(i, 10) == 0
       fprintf('.');
   end
end
fprintf(' Done!\n');

% Close the video file
close(videoObj);

fprintf('\nMovie saved as: %s\n', output_filename);
fprintf('Frame rate: %d fps\n', target_frame_rate);
fprintf('Duration: %.2f seconds\n', output_duration_sec);
fprintf('\nThe MP4 file should play at the correct speed in any video player.\n');