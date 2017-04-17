%Geri Madanguit
%CSE400: Image and Video Processing

% READ ME: To go through a segment of the code, you must change the boolean
% values of 'runAvg', 'runGauss2d', 'runGauss1d', 'runMedian', and
% 'smootherImage'. For the first four values, when one of them is set to 
% true, the rest has to be false. This code is segmented into conditional
% blocks for each task in this given assignment. 
%'runAvg': (problem a) applies averaging by a mean filter. 
%'runGauss2d': (problem b) applies the INT_GAUSS_KER algorithm in a 
% 2D convolution with a 2D gaussian mask. 
%'runGauss1d': (problem b) applies to the SEPAR_FILTER algorithm in two 1D
% covolutions with a 1D gaussian filter. 
%'runMedian': (problem c) applies the median filtering method. 

% This code also includes comments of the
% observation and results of each method and also commented hard code of 
% the third problem, median filtering, in a 3x3 scope- an example of a
% failed version.

clear all;
clc;

noisy1 = imread('NoisyImage1.jpg');
noisy2 = imread('NoisyImage2.jpg');

A=noisy2;

runAvg = false;  %set to TRUE to run AVERAGING BY SMOOTHING FILTER, 
                        % set everything else to FALSE
runGauss2d = false; %set to TRUE to run INT_GAUSS_KER, 
                        % set everything else to FALSE
runGauss1d = false; %set to TRUE to run SEPAR_FILTER, 
                        % set everything else to FALSE
runMedian = true; %set to TRUE to run MEDIAN FILTERING, 
                        % set everything else to FALSE 

smootherImage = true; %set to TRUE to test with NoisyImage1
                       %set to FALSE to test with NoisyImage2

if(smootherImage)
    figure, imshow(noisy1);
    A=noisy1;
else
    figure, imshow(noisy2);
    A=noisy2;
end
    
%A. Apply average smoothing by a mean filter 
% by using two different kernel sizes. Then, compare
% and evaluate your results

if(runAvg)
    for k=1:10
        kernel1 = ones(k);
        smooth = conv2(A,(kernel1)/(k^2),'full');
        smooth = uint8(smooth);
        figure, imshow(smooth);
    end
end

%A. OBSERVATIONS/RESULTS: When using average by smoothing a mean filter,
% 3x3 filter is the 'best' for Noisy1
% 5x5 filter is the 'best' for Noisy2
% the bigger the filter, the more blurrier it gets losing defined lines in
% the image. The edges also become darker, giving a black frame, as the
% kernel size increases

%B. Apply Gaussian smoothing by using two different sigma values
%compare and evaluate your results

%INT_GAUSS_KER

if(runGauss2d)
    
    kerG2r = 3; %gaussian mask row
    kerG2c = 3; %gaussian mask column
    std = 1; %sigma
    
    centerx = round((kerG2r+1)/2);
    centery = round((kerG2c+1)/2);
    
    int_gauss_ker = ones(kerG2r,kerG2c); 
    float_ker = zeros(kerG2r,kerG2c);

    for k=1:kerG2c
        for h=1:kerG2r
            float_ker(h,k) = exp((-(((h-centerx)^2+(k-centery)^2)/2*(std)^2)));
        end
    end
    
    gmin = float_ker(centerx,centery);
    norm_factor = 1/gmin;
    int_gauss_ker = round(norm_factor*float_ker);
    gauss_filtered_im = conv2(int_gauss_ker,A,'full');

    max_gauss_filtered = max(max(gauss_filtered_im));

    gauss_filtered_im = uint8((gauss_filtered_im/max_gauss_filtered)*255);

    figure, imshow(gauss_filtered_im);
end

%B1. OBSERVATIONS/RESULTS: When the row column, and sigma value are all
%equal there was a lack of effect on the image. The larger the gaussian
%mask, the more there is a black frame effect. 
%Noisy1 best = 3x3, sigma =1
%Noisy2 best = 5x5, sigma =1

%SEPAR_FILTER
if(runGauss1d)
    sepC = 5;
    std = 1;
    
    horiz_gauss_mask = ones(1,sepC);
    ctry = round((sepC+1)/2);
   
 
    for col=1:5
        horiz_gauss_mask(1,col) = exp((-(((col-ctry)^2)/2*(std)^2)));
    end

    vert_gauss_mask = horiz_gauss_mask';

    [rows, cols] = size(A);

    for r=1:rows
        ir(r,:) = conv(horiz_gauss_mask, A(r,:));
    end

    [rIR,cIR] = size(ir);

    for c=1:cIR
        res(:,c) = conv(vert_gauss_mask, ir(:,c));
    end

    max_res = max(max(res));
    separ_im = uint8((res/max_res)*255);
    figure, imshow(separ_im);
end

%B2. OBSERVATIONS/RESULTS: 
%Noisy1 best = 1x5, sigma =0.8
%Noisy2 best = 1x3, sigma =1

%C. Apply median filtering by using two different kernel sizes
%{
for j=1:rM
    for k=1:cM
        %leftEdge = noisy1(j,k-1);
        %rightEdge = noisy1(j,k+1);
        %bottomLeft = noisy1(j+1,k-1);
        %bottomEdge = noisy1(j+1,k);
        %bottomRight = noisy1(j+1,k+1);
        middle = noisy2(j,k);
        %topLeft = noisy1(j-1,k-1);
        %topEdge = noisy1(j-1,k);
        %topRight = noisy1(j-1,k+1);
        if j == 1 %toprow 
            if k == 1 %topleft
                A = [middle, noisy2(j,k+1), noisy2(j+1,k), noisy2(j+1,k+1)];
            elseif k == cM %topright
                A = [middle, noisy2(j,k-1), noisy2(j+1,k-1), noisy2(j+1,k)];
            else %topedge
                A = [noisy2(j,k-1), middle, noisy2(j,k+1), noisy2(j+1,k-1), noisy2(j+1,k), noisy2(j+1,k+1)];
            end
        elseif j == rM %bottomrow
            if k == 1 %bottomleft
                A = [middle, noisy2(j-1,k), noisy2(j-1,k+1), noisy2(j,k+1)];
            elseif k == cM %bottomright
                A = [noisy2(j,k-1), middle, noisy2(j-1,k), noisy2(j-1,k-1)];
            else %bottomedge
                A = [noisy2(j,k+1), middle, noisy2(j,k-1), noisy2(j-1,k), noisy2(j-1,k-1), noisy2(j-1,k+1)];
            end
        elseif k == 1 %leftedge
            A = [noisy2(j-1,k), noisy2(j-1,k+1), middle, noisy2(j,k+1), noisy2(j+1,k), noisy2(j+1,k+1)];
        elseif k == cM %rightedge
            A = [noisy2(j+1,k), noisy2(j+1,k-1), middle, noisy2(j,k-1), noisy2(j-1,k-1), noisy2(j-1,k)];
        else
            A = [noisy2(j-1,k-1), noisy2(j-1,k), noisy2(j-1,k+1), noisy2(j,k-1), middle, noisy2(j,k+1), noisy2(j+1,k-1), noisy2(j+1,k), noisy2(j+1,k+1)];
        end
        
        noisy2(j,k) = median(A);  
    end
end
%}      
            %{^^^^^MEDIAN FILTERING FOR 3X3 HARDCODE (correct output, but failed ver.)}%


if(runMedian)
    mfR = 4; %median kernel row
    mfC = 4; %median kernel column
    
    column_to_median = im2col(A, [mfR mfC]);
    median_to_im = median(column_to_median);
    med_filtered_im = col2im(median_to_im,[mfR mfC],size(A));
    figure, imshow(med_filtered_im);
end
%C. OBSERVATIONS/RESULTS:
%This method is very helpful for salt and pepper noise but does little to
%no difference for Noisy1 without destroying the detail.
%Noisy1 best = 3x3
%           2x2 no difference, 4x4 too smooth, too blurred
%Noisy2 best = 4x4

