function fpt = FaceTrack(raw_im, fitting_model, flag)
% use Chehra to locate the 49 facial points for a video

raw_im = im2double(raw_im);
load(fitting_model);

% face detection based on CascadeObjectDetector
faceDetector = vision.CascadeObjectDetector();
faceDetector.ScaleFactor = 2;
bbox = step(faceDetector, raw_im);

if ~isempty(bbox)
    test_init_shape = InitShape(bbox,refShape);
    test_init_shape = reshape(test_init_shape, 49, 2);
    MaxIter=6;
    fpt = Fitting(raw_im,test_init_shape,RegMat,MaxIter);
    
    % show landmarks on facial image
    if flag == 1 && ~isempty(fpt)
        imshow(raw_im, []);
        hold on;
        plot(fpt(:,1),fpt(:,2), 'g*', 'MarkerSize',6);hold off;
    end
else
    fpt = [];
end