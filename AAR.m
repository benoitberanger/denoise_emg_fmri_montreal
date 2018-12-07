function OUT = AAR(Signal, WindowLength, nAvg, verbose)
% Signal       : EMG, EEG, ... with MRI gradient artifacts
% WindowLength : number of data points per TR, sliding window length
% nAvg         : number of averages of the windows around the current window

if ~exist('verbose','var')
    verbose = 1;
end

upscale_factor = 10;


%% Upsample the signal

Signal = interp(Signal,upscale_factor);
WindowLength = WindowLength*upscale_factor;

nWindowInSignal = floor(length(Signal)/WindowLength); % nTR
indexWindow     = 1 : WindowLength : length(Signal);  % iTR

% Memory pre-allocation
carpet_raw   = zeros(nWindowInSignal,WindowLength);
carpet_shift = zeros(nWindowInSignal,WindowLength);
carpet_clean = zeros(nWindowInSignal,WindowLength);


%% Slice the Signal into the WindowLength

for iWindow = 1 : nWindowInSignal
    carpet_raw(iWindow,:) = Signal( indexWindow(iWindow) : indexWindow(iWindow)+WindowLength-1 );
end


%% Compute Timeshift

% Construct Gaussian Kernel for smoothing
% https://fr.mathworks.com/matlabcentral/answers/81689-how-to-implement-convolution-instead-of-the-built-in-imfilter
m=upscale_factor*2; 
n=upscale_factor*2;
sigma=upscale_factor;
[h1,h2]=meshgrid(-(m-1)/2:(m-1)/2, -(n-1)/2:(n-1)/2);
hg= exp(-(h1.^2+h2.^2)/(2*sigma^2));            %Gaussian function
h=hg ./sum(hg(:));

% Smooth
carpet_smoothed = conv2( carpet_raw, h );

% Fetch minimum of the average window

% Select a sub part of the smoothed carpet
lim_1 = WindowLength/2; % middle of the carpet
lim_2 = lim_1 + WindowLength/20; % some points after...
local_carpet_smoothed = carpet_smoothed(:, lim_1:lim_2 ); % look into a smaller carpet, to avoid the edges 

% Compute the mean, and fetch the minimum : this will be our reference point for the timeshift
mean_local_carpet_smoothed = mean(local_carpet_smoothed,1);
[~,I]= min(mean_local_carpet_smoothed); % Where should I compute the timeshift ?
I = I + lim_1; % real index of the minimum from the whole carpet

% if verbose > 0
%     figure('Name','local Carpet smoothed, used to fetch where the timeshift will be computed')
%     ax(1) = subplot(4,1,1:3);
%     image(local_carpet_smoothed,'CDataMapping','scaled')
%     colormap(gray(256))
%     colorbar
%     ax(2) = subplot(4,1,4);
%     image(mean_local_carpet_smoothed,'CDataMapping','scaled')
%     colormap(gray(256))
%     colorbar
%     linkaxes(ax,'x')
% end

% Reduce our carpet around this supposed local minimum
local_carpet = carpet_smoothed(:,I-10*upscale_factor/2:I+10*upscale_factor)/2;
if verbose > 0
    figure('Name','Carpet around the minimum, where timeshift is computed')
    subplot(2,1,1)
    image(local_carpet,'CDataMapping','scaled')
    colormap(gray(256))
    colorbar
end

% Fetch the local minimum for each window
local_min = zeros(nWindowInSignal,1);
for iWindow = 1 : nWindowInSignal
    [ ~, local_min(iWindow) ] =  min(local_carpet(iWindow,:));
end

% This the number sample to shift to realign each window
shift_value = local_min(1) - local_min;
if verbose > 0
    subplot(2,1,2)
    plot(shift_value,'LineWidth',2);
end


%% Is it necessary to do a timeshift ?

x = 0:length(shift_value)-1;
y = shift_value';
p = polyfit(x,y,1);
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
if verbose > 0
    hold on
    title(sprintf('R2 = %g', rsq))
    plot(x,yfit)
end


%%  Apply shift (or not)

if rsq>0.9
    
    for iWindow = 1 : nWindowInSignal
        if iWindow == 1
            prev = carpet_raw(iWindow,:);
        else
            prev = carpet_raw(iWindow-1,:);
        end
        curr     = carpet_raw(iWindow  ,:);
        concat                  = [prev curr];
        concat_shifted          = circshift( concat , [0 shift_value(iWindow)] );
        carpet_shift(iWindow,:) = concat_shifted(WindowLength+1:WindowLength*2);
    end
    
else
    
    carpet_shift = carpet_raw;
    
end


%% AAR


for iWindow = 1 : nWindowInSignal
    
    carpetSelection = zeros(nAvg,size(carpet_shift,2));
    counter = 0; % line counter
    for iSlidingWindow = iWindow-(nAvg-1)/2 : iWindow+(nAvg-1)/2
        if iSlidingWindow <= 0
            % pass
        elseif iSlidingWindow >= nWindowInSignal
            % pass
        else
            counter = counter + 1;
            carpetSelection(counter,:) = carpet_shift(iSlidingWindow,:);
        end
    end
    
    meanCarpet = sum(carpetSelection,1)/counter;

    carpet_clean(iWindow,:) = carpet_shift(iWindow,:) - meanCarpet;
    
end

if verbose > 0
    
    figure('Name','Last TRs and their avg')
    hold on
    plot(meanCarpet,'-r','LineWidth',4,'DisplayName','avg')
    plot(carpetSelection(1:counter,:)')
    legend show
    
    figure('Name','Carpet plot')
    
    ax(1) = subplot(3,1,1);
    image(carpet_raw,'CDataMapping','scaled')
    colormap(gray(256))
    colorbar
    title('raw')
    
    ax(2) = subplot(3,1,2);
    image(carpet_shift,'CDataMapping','scaled')
    colormap(gray(256))
    colorbar
    title('after timeshift')
    
    ax(3) = subplot(3,1,3);
    image(carpet_clean,'CDataMapping','scaled')
    colormap(gray(256))
    colorbar
    title('final : after AAR')
    
    linkaxes(ax,'xy')
    
end


%% Prepare output

OUT = carpet_clean';
OUT = OUT(:);                          % change from carpet signal (1 vect)
OUT = [OUT' Signal(numel(OUT)+1:end)]; % append the last samples not treated
OUT = downsample(OUT,upscale_factor);


end % function
