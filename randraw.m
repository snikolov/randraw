% Load an image to randomly draw
name = 'circle.png';

CREATE = 1;
IMAGE_FIGURE = 1;

if IMAGE_FIGURE
  figure
end

if CREATE
  I = imread(name);

  % Convert to grayscale
  if numel(size(I)) == 3
    I = rgb2gray(I);
  end

  I = imresize(I, 0.4);
  % Create energy function
  th = 130;
  Ith = double(I < th);

  if IMAGE_FIGURE
    imagesc(Ith);
    colormap gray
    pause;
  end

  Iblr = conv2(conv2(Ith, ones(1, 25), 'same'), ones(25, 1), 'same');

  Iphi = mask2phi(Ith);
  Iphi = -Iphi;
  Iphi = Iphi - min(Iphi(:));

  % Add blur component
  I = Iphi + 0.0 * mean(Iphi(:)) * Iblr / mean(Iblr(:));
  I = I - min(I(:)) + 0.1;
  
  save([name, '.mat'], 'I');
else
  % Load energy function
  load([name, '.mat']);
end

if IMAGE_FIGURE
  surf(double(imresize(I, 0.5)))
  colormap gray
  pause;
end

% Brush parameters
brushRadius = 4;
sigmaBrush = 0.2 * brushRadius;
[x y] = meshgrid(-brushRadius:brushRadius);
brush = exp(-(x.^2 + y.^2)/(2 * sigmaBrush^2));

% Initialize N random walkers
rand('twister',6657)
N = 10;
[r, c] = size(I);
rwr = randi(r, N, 1);
rwc = randi(c, N, 1);
% Proposal distribution parameters.
sigmaQ = sigmaBrush * 1;

% Time steps.
T = 100000;

% Drawing buffer.
D = zeros(size(I));

DRAWING_FIGURE = 1;
if DRAWING_FIGURE
  figure
  subplot(211)
  imagesc(I)
  hold on;
end

valSums = [];
plotMod = 5;
for t = 1:T
  % Draw
  for ridx = 1:numel(rwr)
    if rwr(ridx) > brushRadius && rwr(ridx) <= r - brushRadius ...
	  && rwc(ridx) > brushRadius && rwc(ridx) <= c - ...
	  brushRadius
      D(rwr(ridx) - brushRadius:rwr(ridx) + brushRadius, rwc(ridx) ...
	- brushRadius:rwc(ridx) + brushRadius) = brush + D(rwr(ridx) - brushRadius:rwr(ridx) + brushRadius, rwc(ridx) ...
	- brushRadius:rwc(ridx) + brushRadius);
    end
  end
  if DRAWING_FIGURE
    subplot(211)
    scatter(rwc, rwr);
    %pause;
    %imagesc(-D)
    colormap gray;
    if ~mod(t, plotMod)
      subplot(212)
      valSums(end + 1) = sum(diag(I(rwr, rwc)));
      plot(valSums)
      drawnow;
    end
    title(sprintf('t=%d', t));
  end
  % Proposed random step
  prwr = mod(round(rwr + normrnd(0, sigmaQ, N, 1)) - 1, r) + 1;
  prwc = mod(round(rwc + normrnd(0, sigmaQ, N, 1)) - 1, c) + 1;

  for pidx = 1:numel(rwr)
    new2old = I(prwr(pidx), prwc(pidx)) / I(rwr(pidx), rwc(pidx));
    accepted = new2old >= 1 || rand < new2old;
    %scatter(prwc(pidx), prwr(pidx), 12, 'square', 'filled');
    %fprintf('new2old: %.4f, accepted: %d\n', new2old, accepted);
    
    if accepted
      if norm([rwc(pidx), rwr(pidx)] - [prwc(pidx), prwr(pidx)]) < ...
	    6 * sigmaQ
	%line([rwc(pidx), prwc(pidx)], [rwr(pidx), prwr(pidx)]);
      end
      rwr(pidx) = prwr(pidx);
      rwc(pidx) = prwc(pidx);
    end
  end
  rwc = prwc;
  rwr = prwr;
end
