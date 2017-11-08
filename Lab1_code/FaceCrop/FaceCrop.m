function norm_face = FaceCrop(raw_im, ptdat, normal)
% face crop according to the landmark
% In this version, we estimate the left eye coordinate according to the No. 21,
% 22, 24 and 25 landmarks provided by Chehra tracker. Right eye coordinate
% according to No. 27, 28, 30 and 31.

if size(raw_im, 3) == 3
    raw_im = rgb2gray(raw_im);
end

% Extract the eye coordinate
p1 = ptdat(21, :);
p2 = ptdat(22, :);
p3 = ptdat(24, :);
p4 = ptdat(25, :);

leye = (p1 + p2 + p3 + p4)./4; % the coordinate of left eye

p1 = ptdat(27, :);
p2 = ptdat(28, :);
p3 = ptdat(30, :);
p4 = ptdat(31, :);

reye = (p1 + p2 + p3 + p4)./4; % the coordinate of right eye

x1 = leye(1, 1);
y1 = leye(1, 2);
x2 = reye(1, 1);
y2 = reye(1, 2);

% distance between left and right eyes
DisBetEye = sqrt ((y1 - y2) * (y1 - y2)+ (x1 - x2)*(x1 - x2));
sina= (y1-y2)/DisBetEye; % the rotation angle
cosa = (x2-x1)/DisBetEye;
lefttopy = y1 + DisBetEye * 0.4 * sina -DisBetEye * 0.6 * cosa;
lefttopx = x1 - DisBetEye * 0.4 * cosa - DisBetEye * 0.6 * sina;

% pre-define the size of cropped face image
faceHeight = round(DisBetEye * 2.2);
faceWidth = round(DisBetEye * 1.8);
norm_face = zeros(faceHeight, faceWidth);
[wi, hi] = size(raw_im);

for h = 1:1:faceHeight
    starty = lefttopy + h * cosa;
    startx = lefttopx + h * sina;
    for w = 1:1:faceWidth
        if uint16(starty - w * sina) > wi
            norm_face(h,w) = raw_im(uint16(wi), uint16(startx + w * cosa));
        elseif uint16(startx + w * cosa) > hi
            norm_face(h,w) = raw_im(uint16(starty - w * sina), uint16(hi));
        else
            norm_face(h,w) = raw_im(uint16(starty - w * sina), uint16(startx + w * cosa));
        end
    end
end

if normal == 1
    norm_face = imresize(uint8(norm_face), [128 128]); % keep the same distance of eyes for all images
end




