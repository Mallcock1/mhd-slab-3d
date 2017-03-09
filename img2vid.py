
import os

# jsallcock 11/16

def image2video(filepath=None, prefix='', in_extension='png',
            out_extension='avi', output_name=None, fps=10, n_loops=1, 
            delete_images=False, delete_old_videos=False, res=1080,
            overlay=True, cover_page=False):
                
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
    
    delete_old_videos_cmd = 'DEL "'+filepath+output_name+'.'+out_extension+'" \
    "'+filepath+output_name+'_overlay.'+out_extension+'" "'+filepath+output_name+'_overlay2.'+out_extension+'"'


    logo_height = str(int(100. * res / 1080.))
    logo_width = str(int(300. * res / 1080.))
    
    overlay_image_sp2rc = 'ffmpeg -y -i '+filepath+output_name+'.'+out_extension+' -i \
    sp2rc_logo2.png -filter_complex "[1:0] scale=-1:'+logo_height+' [logo]; \
    [0:0][logo] overlay=main_w-'+logo_width+':main_h-'+logo_height+'" -c:a copy \
    '+filepath+output_name+'_overlay.'+out_extension

    overlay_image_swat = 'ffmpeg -y -i '+filepath+output_name+'_overlay.'+out_extension+' -i \
    swat_logo2.png -filter_complex "[1:0] scale=-1:'+logo_height+' [logo]; \
    [0:0][logo] overlay=main_w-overlay_w:main_h-'+logo_height+'" -c:a copy \
    '+filepath+output_name+'_overlay2.'+out_extension        



# create .txt file to loop through
    loop_list = "(for /l %i in (1,1,"+str(n_loops)+") do @echo file \
    '"+filepath+output_name+'_overlay2.'+out_extension+"') > "+filepath+'mylist.txt'
    
# loop vid n_loops times
    loop = 'ffmpeg -y -f concat -safe 0 -i '+filepath+'mylist.txt -c copy \
    '+filepath+output_name+'_overlay2_looped.'+out_extension

    
    os.system(image_2_video)
    if overlay == True:
        os.system(overlay_image_sp2rc)
        os.system(overlay_image_swat)
    if n_loops!=1:
        os.system(loop_list)
        os.system(loop)
    if delete_images == True:
        os.system(delete_images_cmd)
    if delete_old_videos == True:
        os.system(delete_old_videos_cmd)
    

#img2vid(prefix=prefix, output_name='video', out_extension='mp4', fps=20, n_loops=4, delete_images=True)