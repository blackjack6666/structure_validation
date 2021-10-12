import os
import moviepy.video.io.ImageSequenceClip

fps=1

image_folder = 'D:/data/alphafold_pdb/native_digest_time_laps'
image_files = [image_folder+'/'+img for img in os.listdir(image_folder) if 'HSPE1' in img]
print (image_files)
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('D:/data/alphafold_pdb/native_digest_time_laps/HSPE1.mp4')