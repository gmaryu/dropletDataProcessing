function averageIntensityPerLabel = getAverageIntensityPerLabel(image, labeledImage)
    %    getAverageIntensityPerLabel averages the pixels in the image using the labels in labeledImage
    %
    %    averageIntensityPerLabel = getAverageIntensityPerLabel(image, labeledImage)
    %
    %    Input:
    %        image (MxN) double - The image to be averaged
    %        labeledImage (MxN) numeric - The image with labels
    %
    %    Output:
    %        averageIntensityPerLabel (2xL) double - Two rows, one for the label, one for the average intensity
    arguments
        image (:,:) double
        labeledImage (:,:) {mustBeNumeric}
    end
    props = regionprops(labeledImage, image, 'MeanIntensity');
    averageIntensityPerLabel = [cat(1, props.MeanIntensity); cat(1, props.Label)]';
end