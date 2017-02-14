
import os
#from matplotlib import pyplot as plt
#from matplotlib import image as mpimg
#import matplotlib.colors
#import glob

# jsallcock 11/16

def image2video(filepath=None, prefix='', in_extension='png',
            out_extension='avi', output_name=None, fps=10, n_loops=1, 
            delete_images=False, delete_old_videos=False, res=1080):
    if filepath == None:
        filepath = 'D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_animations\\'
#        filepath = u'/media/matthew/W7_DATA/my_work/projects/Asymmetric_slab/Python/visualisations/ffmpeg/'
    if output_name == None:
        output_name = prefix
    in_fps = fps
    out_fps = fps
    start_frame = 1

    image_2_video = 'ffmpeg -y -framerate '+str(in_fps)+' -start_number '+str(start_frame)+' -i \
    '+filepath+prefix+'%d.'+in_extension+' \
    -c:v libx264 -r '+str(out_fps)+' -pix_fmt yuv420p \
    '+filepath+output_name+'.'+out_extension

    delete_images_cmd = 'DEL "'+filepath+prefix+'*.'+in_extension+'"'
    
    delete_old_videos_cmd = 'DEL "'+filepath+prefix+'.'+out_extension+'" \
    "'+filepath+prefix+'_overlay.'+out_extension+'" "'+filepath+prefix+'_overlay2.'+out_extension+'"'

#    image_2_video = 'ffmpeg -framerate '+str(in_fps)+' -s 1920x1080 -start_number '+str(start_frame)+' -i \
#    '+filepath+prefix+'.'+in_extension+' \
#    -r '+str(out_fps)+' -vcodec libx264 -crf 25 -pix_fmt yuv420p \
#    '+filepath+output_name+'.'+out_extension

##THIS WORKS< OVERLAYING ONE LOGO
#    overlay_image_sp2rc = 'ffmpeg -i '+filepath+output_name+'.'+out_extension+' -i \
#    sp2rc_logo2.png -filter_complex "[1:0] scale=-1:100 [logo]; \
#    [0:0][logo] overlay=main_w-300:main_h-100" \
#    -pix_fmt yuv420p -c:a copy \
#    '+filepath+output_name+'_overlay.'+out_extension
#
#    overlay_image_swat = 'ffmpeg -i '+filepath+output_name+'_overlay.'+out_extension+' -i \
#    swat_logo2.png -filter_complex "[1:0] scale=-1:100 [logo]; \
#    [0:0][logo] overlay=main_w-overlay_w:main_h-100" \
#    -pix_fmt yuv420p -c:a copy \
#    '+filepath+output_name+'_overlay2.'+out_extension    
    logo_height = str(int(100. * res / 1080.))
    logo_width = str(int(300. * res / 1080.))
    
#THIS WORKS< OVERLAYING ONE LOGO
    overlay_image_sp2rc = 'ffmpeg -y -i '+filepath+output_name+'.'+out_extension+' -i \
    sp2rc_logo2.png -filter_complex "[1:0] scale=-1:'+logo_height+' [logo]; \
    [0:0][logo] overlay=main_w-'+logo_width+':main_h-'+logo_height+'" -c:a copy \
    '+filepath+output_name+'_overlay.'+out_extension

    overlay_image_swat = 'ffmpeg -y -i '+filepath+output_name+'_overlay.'+out_extension+' -i \
    swat_logo2.png -filter_complex "[1:0] scale=-1:'+logo_height+' [logo]; \
    [0:0][logo] overlay=main_w-overlay_w:main_h-'+logo_height+'" -c:a copy \
    '+filepath+output_name+'_overlay2.'+out_extension        

#    overlay_image = 'ffmpeg -i '+filepath+output_name+'.'+out_extension+' -i \
#    sp2rc_logo.png -i swat_logo.png -filter_complex "[1:0] scale=100:100 [sp2rc-logo]; \
#    [2:0] scale=100:100 [swat-logo]; [0:0][sp2rc-logo] overlay=main_w-100:main_h-100" \
#    -pix_fmt yuv420p -c:a copy \
#    '+filepath+output_name+'_overlay.'+out_extension


#    [1:v] scale=20:20, 
# create .txt file to loop through
    loop_list = "(for /l %i in (1,1,"+str(n_loops)+") do @echo file \
    '"+filepath+output_name+'_overlay2.'+out_extension+"') > "+filepath+'mylist.txt'
    
# loop vid n_loops times
    loop = 'ffmpeg -y -f concat -safe 0 -i '+filepath+'mylist.txt -c copy \
    '+filepath+output_name+'_overlay2_looped.'+out_extension

#Below doesnt work
#    loop = 'ffmpeg -stream_loop '+str(n_loops)+' -i '+filepath+output_name+'.'+out_extension+'\
#    -c copy '+filepath+output_name+'.'+out_extension
    
    os.system(image_2_video)
    if delete_images == True:
        os.system(delete_images_cmd)
    os.system(overlay_image_sp2rc)
    os.system(overlay_image_swat)
    if n_loops!=1:
        os.system(loop_list)
        os.system(loop)
    

#img2vid(prefix=prefix, output_name='video', out_extension='mp4', fps=20, n_loops=4, delete_images=True)