# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/make_colortable.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import sys
import os
try:
    import ImageFont, ImageDraw, Image
except ImportError:
    from PIL import ImageFont, ImageDraw, Image

__general = r'/home/arthur/globalnutrients/generalcode/trunk'
if os.path.exists(__general): 
    sys.path.insert(0, __general)
    print(__general + " is added to the python search path for modules.") 

import ascraster
import my_sys

# Give name of the output file
base = "colortable"

# Create a grid of 16 by 16
grid = ascraster.Asciigrid(nrows=16,ncols=16, xllcorner = 0, yllcorner = 0, cellsize = 1)
# Fill data of grid with values (0.5 upto 255.5)
for icell in range(grid.length):
    grid.set_data(icell,icell+0.5)
maxV = max(grid.values)
grid1 = grid.resize(75)
grid1.save_as_bitmap(base + ".jpg",minV=0.0,maxV=maxV,mode="RGB",color=(1,1,1),number_of_classes = 255, classtype="colorset_general.txt",\
                    nodatacolor=(225,225,225))
# The rest of the file is for purpose of making videos. Not used here
sys.exit(0)

def slider(im, left, length, pos, min_pos, max_pos, pos_text = False):
    '''Adds a time slider with current position to a image.
    
    @param im: Image instance to draw the slider on
    @type im: Image instance
    @param left: Coordinate left side of slider line
    @type left: tuple/list
    @param length: Length of slider line
    @type length: integer
    @param pos: position of frame
    @type pos: integer/float
    @param min_pos: Start position
    @type min_pos: integer/float
    @param max_pos: End position
    @type max_pos: integer/float
    @param pos_text: show pos text above current position
    @type pos_text: boolean
    Example: slider(img, (200,1070), 1800, 5*(num-1)+1900, 1900, 2000, pos_text = True)
    '''
    draw = ImageDraw.Draw(im)
    draw.line([left, (left[0]+length, left[1])], fill=(0,110,154), width=4)

    center = (length * (pos - min_pos)) / (max_pos - min_pos) + left[0]
    radius = 10
    bbox = (center-radius,left[1]-radius,center+radius,left[1]+radius)
    draw.ellipse(bbox, fill=(157,0,110))
    if pos_text:
        font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSansCondensed-BoldOblique.ttf", 50)
        draw.text((center-52,left[1]-55),str(pos),(0,110,154),font=font)
    del draw

def do_it(listelem):
        filename = listelem[0]
        num = listelem[1]
        base = "lala"
        if (num < 10):
            base = base + "000" + str(num)
        elif (num < 100):    
            base = base + "00" + str(num)
        elif (num < 1000):    
            base = base + "0" + str(num)
        else:
            base = base + str(num)

        grid1 = ascraster.Asciigrid(ascii_file=filename)

        isogrid = ascraster.Asciigrid(ascii_file="../gcountry.map")
        for icell in range(isogrid.length):
            iso = isogrid.get_data(icell)
            if (iso != None):
                # Make Greenland nodata
                if (int(float(iso)) == 304): 
                    grid1.set_data(icell,0.0)

        grid = grid1.resize(3)
        #grid = grid1
        #grid.write_ascii_file(str(num) +".asc")

        #grid.save_as_bitmap("qq_"+str(num)+".png",minV=0.0001,maxV=maxV,mode="RGB",color=(1,0,0),number_of_classes = 150, classtype="log")
        grid.save_as_bitmap(base + ".png",minV=0.0,maxV=maxV,mode="RGB",color=(1,1,1),number_of_classes = 100, classtype="colorset.txt",\
                            nodatacolor=(225,225,225))

        # Put year in the image

        img = Image.open(base+".png")
        draw = ImageDraw.Draw(img)
        font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSansCondensed-BoldOblique.ttf", 60)
        draw.text((10, 30), str(5*(num-1)+1900), (0,0,0),font=font)
        font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSansCondensed-BoldOblique.ttf", 30)
        draw.text((10, 980), "pbl.nl", (100 , 100,  100),font=font)

        slider(img, (200,1070), 1800, 5*(num-1)+1900, 1900, 2000, pos_text = True)

        img.save(base+".png")


        print(filename + " is ready.")
        #os.system('convert qq_'+str(num)+'.png -size 360X720 -scale 200% ' + base + '.png')
        #my_sys.my_removefile('qq_'+str(num)+'.png')

# Find maximum scale for all the maps.

maxnum = -1
#Make a loop over the years.
for iyear in range(1850,2100):
    if (os.path.isfile(os.path.join(outputdir,os.path.join(str(iyear),basename)))):
        # Make base number so that it is fixed format
        maxfile = iyear

# Get maximum value for the last map
grid = ascraster.Asciigrid(ascii_file=os.path.join(outputdir,os.path.join(str(maxfile),basename)))
maxV = max(grid.values)

list1=[]
num = 1
for iyear in range(1850,2100):
    if (os.path.isfile(os.path.join(outputdir,os.path.join(str(iyear),basename)))):
        list1.append([os.path.join(outputdir,os.path.join(str(iyear),basename)),num])
        num += 1

p = multiprocessing.Pool()
p.map(do_it,list1)

os.system('convert -delay 60  *.png qq.mp4')
# Set bitrate
os.system('ffmpeg -i qq.mp4 -y -b 1200 filmpje.mpg')
    


