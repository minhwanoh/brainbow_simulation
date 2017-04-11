function writeSimulatedImageSeq()
tempVol = overallRawVolume;
data = uint32(overallRawVolume);

for i=1:100;
    outputFileName = ['simulated_data/simVol_', num2str(i), '.tif'];;
    % This is a direct interface to libtiff
    t = Tiff(outputFileName,'w');

    % Setup tags
    % Lots of info here:
    % http://www.mathworks.com/help/matlab/ref/tiffclass.html
    tagstruct.ImageLength     = size(data,1);
    tagstruct.ImageWidth      = size(data,2);
    tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample   = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip    = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software        = 'MATLAB';
    t.setTag(tagstruct)
    
    t.write(data(:,:,i));
    t.close();
end