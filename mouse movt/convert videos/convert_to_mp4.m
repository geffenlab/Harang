function convert_to_mp4( reader, filepath )
%convert_to_mp4 convert video in reader (VideoReader) to mp4 @ filepath

v = VideoWriter(filepath, 'MPEG-4');

open(v)
while hasFrame(reader)
    writeVideo(v, readFrame(reader));
end
close(v)

end

