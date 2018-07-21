clc;

%reads in image and converts to greyscale
O = imread('Zebra.jpg');
G = rgb2gray(O);

%Variable declarations
[x,y] = size(G); %size of original image
X = 3*x; %scale of 3 so new X/Y is *3 larger
Y = 3*y;
RX = (x/X);%ratio small to large image
RY = (y/Y);


%declare output size
nnOut = zeros(X,Y, class(G));
blOut = zeros(X,Y, class(G));

%nested loop to iterate through output image pixels
for i=1:X
    %stores a float value for x for pixel
    ii = (RX * i) + (0.5 * (1 - 1/3));
    for j=1:Y
        %stores a float value for y for pixel
        jj = (RY * j) + (0.5 * (1 - 1/3));
        
        %getting x1, x2
        jj(jj < 1) = 1; %error handling
        jj(jj > y - 0.001) = y - 0.001; %ensures always in range
        x1 = floor(jj);
        x2 = x1+1;
        
        %Getting y1, y2
        ii(ii < 1) = 1; 
        ii(ii > x - 0.001) = x - 0.001;
        y1 = floor(ii);
        y2 = y1+1;
        
        %stores intensity values at 4 neighbours
        N11 = G(y1,x1);
        N21 = G(y2,x1); 
        N12 = G(y1,x2);
        N22 = G(y2,x2);
        
        %linear interpolations in the x direction
        R1 = ((x2 - jj)/(x2 - x1))*N11 + ((jj - x1)/(x2 - x1))*N12;
        R2 = ((x2 - jj)/(x2 - x1))*N21 + ((jj - x1)/(x2 - x1))*N22;
        
        %linear interpolation in the y direction and combines for bilinear
        blOut(i,j) = ((y2 - ii)/(y2 - y1))*R1 + ((ii - y1)/(y2 - y1))*R2;
        
        %gets nearest neighbour to pixel
        ii_NN = round( (i-1)*(x-1)/(X-1)+1 );
        jj_NN = round( (j-1)*(y-1)/(Y-1)+1 );

        %assigns neighbour to appropriate pixel in output
        nnOut(i,j) = G(ii_NN, jj_NN);
        
    end
end

%displays outputs
colormap gray;
subplot(1,3,1), imagesc(G), title("Input image");
subplot(1,3,2), imagesc(blOut), title("Bilinear Interpolation");
subplot(1,3,3), imagesc(nnOut), title("Nearest Neighbour Interpolation");

clear; close all; clc; 

%read in input image
G=imread('SC.png');
%variable declarations
[m,n]=size(G); %row and column length
A=80; 
B=100;
C=220;
%create new matrix for output of class G (input)
O=zeros(m,n,class(G));
%nested loop to iterate through all input matrix elements
for i=1:m
    for j=1:n
        %reassigns pixel values in geyscale range 80 - 100 to 220
        if G(i,j)>=A && G(i,j)<=B
            O(i,j)=C;
        %maps all other pixels to new matrix
        else
            O(i,j)=G(i,j);

        end
    end
end
%displays outputs
colormap gray;
subplot (1,2,1), imagesc(G), title("Original");
subplot(1,2,2), imagesc(O), title("Piecewise-Linear transformation");

clear; close all; clc; 

%read in input image and convert to greyscale
I = imread('Noisy.png');
G = rgb2gray(I);
%adds a 2 pixel padding
O = padarray(G, [2,2],'replicate');
%convert greyscale values to double
X = im2double(O);
%variable declarations
[m ,n] =size(G); %row and column length
SE =strel('square', 5); % 5*5 square structuring element
%declare new matricies for output
newImage= zeros(m,n); %output matrix for mean filter
finalMedian = zeros(m,n); %output matrix for median filter
medianImage= zeros(1,25); %matrix to store neighbourhood of each pixel

%nested loop to iterate through input matrix
for i = 1:m
    for j = 1:n
        %resets sum variable and counter
        sum = 0;
        c = 1;
        %nested loop to find the neighbourhood of a given pixel
        for k = i:i+4  
            for l = j:j+4
                %calculates mean 
                sum = sum + (1/25)*X(k,l);
                %stores neighbourhood of pixel in 1d array
                medianImage(1,c) =  X(k,l);
                %increments counter
                c = c+1;
            end
        end
    %calculates median for neighbourhood and stores to output matrix    
    finalMedian(i,j) = median(medianImage);    
    %stores mean for a pixel in output matrix
    newImage(i,j) = sum;    
    end
    
end    
%displays outputs   
colormap gray;
subplot (1,3,1), imagesc(O), title("Original");
subplot(1,3,2), imagesc(newImage), title("Mean Filter");
subplot(1,3,3), imagesc(finalMedian), title("Median Filter");
