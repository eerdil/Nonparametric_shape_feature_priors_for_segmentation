clear all
close all
clc

testImageName = 'testImage5';
testImageFolderName = 'TestImage';

load(sprintf('%s/%s', testImageFolderName, testImageName)); % load test image

numberOfShapesInEachClass = 10; % number of shape in each shape class in training set
numberOfClasses = 10; % number of shape classes in training set

%% providing initial curve
sz_i = size(testImage, 1);
sz_j = size(testImage, 1);
figure, imagesc(testImage); colormap(gray); axis('off');%display test image to ask initial curve
Psi = initialLevelSet(20, sz_i, sz_j);
hold on,
contour(Psi, [0 0], 'r'); % display initial curve on the figure
hold off;

%% construct level set representation of shapes in training set
load('trainingShapesAndFeatures/shapeTraining.mat');
trainingIMatrix = AlignedShapeMatrix;
numberOfShapesInTrainingSet = size(AlignedShapeMatrix, 2);
trainingPhiMatrix = zeros(sz_i * sz_j, numberOfShapesInTrainingSet);
for i = 1:numberOfShapesInTrainingSet
    curShape = double(reshape(trainingIMatrix(:, i), [sz_i sz_j]) > 0);
    trainingPhiMatrix(:, i) = reshape(FMM(0.5 - curShape, 0), [sz_j * sz_i, 1]);
end

%% curve evolution with data term
numberOfIterations = 200;
pose = zeros(1, 4); pose(4) = 360 / 360;
dt = 1; % gradient step size
alpha = 1; % weight of shape force
beta = 0;

for i = 1:numberOfIterations
    Evolve_data(double(testImage), trainingIMatrix, trainingPhiMatrix, Psi, pose, 5, dt, ...
        alpha, beta, 1, numberOfShapesInEachClass * numberOfClasses);

    DrawContoursAndMassCenters(testImage, Psi, sz_i);
end
hold off;

%% Feature extraction and loading feature priors
load('trainingShapesAndFeatures/featureTraining.mat'); % load feature training set
featureTraining = featureTraining / 255;
temp = Psi < 0; % find foreground of the curve
extractedFeature = mean(testImage(find(temp == 1)))/255; % take the mean intensity inside the curve as feature

%% Segmentation using nonparametric shape and feature priors

numberOfIterations = 200; % number of iterations for curve evolution

dt = 0.00001; % gradient step
beta = 5; % weight of the shape term
alpha = 1; % weight of the data term
narrowBandRadius = 5; % radius of narrow band for the narrow band representation of the curve
display = 1; % 1 - display curve evolution, 0 - do not display

kernelSizeShape = zeros(1);
kernelSizeShape = computeKernelSizeShape(trainingPhiMatrix, numberOfClasses, numberOfShapesInEachClass, kernelSizeShape, sz_i, sz_j);
kernelSizeFeature = zeros(1);
kernelSizeFeature = computeKernelSizeFeature(featureTraining, numberOfClasses, numberOfShapesInEachClass, kernelSizeFeature, 1, 100);

for i = 1:numberOfIterations
    shapeF = zeros(size(testImage));
    Evolve_ver12L2(double(testImage), trainingIMatrix, trainingPhiMatrix, Psi, pose, ...
        narrowBandRadius, dt, alpha, beta, 1, numberOfShapesInEachClass * numberOfClasses, ...
        kernelSizeShape, kernelSizeFeature, featureTraining, extractedFeature, shapeF);
    if(display && mod(i, 1) == 0)
        imagesc(testImage); colormap(gray); axis('off');
        title(sprintf('iteration %d / %d', i, numberOfIterations));
        hold on;
        contour(Psi, 'LineWidth', 3, 'LineColor', [1 0 0], 'LevelList', 0);
        drawnow;
        hold off;
    end
end


















