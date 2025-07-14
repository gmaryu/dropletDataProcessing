function retv = convertAreaPixelsToVolume(input, pixelToLength)
    retv = 4 / 3 * pi * power(sqrt(input / pi) * pixelToLength, 3);
end