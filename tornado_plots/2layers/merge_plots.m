kxx = imread('kxx_1layer.png');
kyy = imread('kyy_1layer.png');
kzz = imread('kzz_1layer.png');

% Create a new figure
figure('Units', 'pixels', 'Position', [400, 350, 3000, 400]);

% Display the first image in the first subplot
subplot(1, 3, 1);
imshow(kxx);

% Display the second image in the second subplot
subplot(1, 3, 2);
imshow(kyy);

% Display the third image in the third subplot
subplot(1, 3, 3);
imshow(kzz);

% Save the figure as an image file
print(gcf, 'merged_images.png', '-dpng', '-r800');

% Close the figure
close(gcf);