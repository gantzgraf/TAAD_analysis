''' A collection of functions which alter or add elements to a matplotlib/seaborn plot. '''
import matplotlib as mpl
# allows one to run matplotlib and seaborn headless
mpl.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import PIL
import textwrap
from PIL import ImageDraw, ImageFont, Image

def line_between_plots(axs, x1, x2, height, string, fontsize=12, extend=2):
    ''' Draw a horizontal line and string between two points with 
        vertical lines that extend at said points e.g.
 
                      p=0.02
                 ---------------
                |               |
 
    Args:
        axs: list of axes to draw upon
        x1: start of horizontal line
        x2: end of horizontal line
        height: where on the y-axis the horizonatl line be drwan
        string: text to place in the middle of the horizontal line
        fontsize: fontsize of the string
        extend: points to extend the vertical lines by

    Notes:
        if this drawing is not visible on your axes then you may
        have to manually set the ylim to ensure it is visible.
    '''
    if not isinstance(axs, list):
        axs = [axs]
        
    for ax in axs:
        ax.add_line(Line2D(xdata=[x1, x2], ydata=[height, height], color='black'))
        ax.add_line(Line2D(xdata=[x1, x1], ydata=[height-extend, height], color='black'))
        ax.add_line(Line2D(xdata=[x2, x2], ydata=[height-extend, height], color='black'))

        # add string in the middle of horizontal line
        between_points = (x2+x1)/2
        ax.text(x=between_points, y=height+(height*0.01), 
                s=string, fontsize=fontsize,
                horizontalalignment='center')
        return ax

def rename_xtick(df, category_name, counts=True, NA=False, validated_path=False):
    ''' iterate over category counts to get the count for 
        each category name and bind to the xtick_label. If the 
        name is in the name2label dict, then change the name.
        
        counts will not be added to the label if the var counts is False
        BE CARFUL: ensure ticks are renamed as expected
    '''
    name2label={
        'unknown': "Unknown", 'yes': "Family History", 'no': "No Family History"
    }
        
    if NA:
        df = df.dropna(subset=[category_name])
    # filter out pathogenic and likely pathogenic variants that have not been validated by sanger sequencing
    if validated_path is True:
        df = df[(((df['Category'] == "Pathogenic") | (df['Category'] == "Likely Pathogenic")) & 
                 (df['validation'] == 1)) | ((df['Category'] != "Pathogenic") & (df['Category'] != "Likely Pathogenic"))]
    # get the counts from the category and concatenate to the category name
    sizes = df.groupby(category_name).size().iteritems()
    # decide whether to add counts to the xtick labels
    if counts is True:
        name_count = [(str(name)+"\n"+"n = "+str(n)) for name, n in sizes]
    else:
        name_count = [(str(name)+"\n") for name, n in sizes]
    # decide whether to replace the name portion of name_count
    name_change = []
    for info in name_count:
        name = info.split("\n")[0]
        if name in list(name2label.keys()):
            new_name = info.replace(name, name2label.get(name))
            name_change.append(new_name)
        else:
            name_change.append(info)
    # Warn user        
    print("\n\nBE CARFUL: ensure xticks are renamed as expected!\n\n")                    
    return name_change

def RoundToSigFigs( x, sigfigs ):
    """
    Rounds the value(s) in x to the number of significant figures in sigfigs.

    Restrictions:
    sigfigs must be an integer type and store a positive value.
    x must be a real value or an array like object containing only real values.
    """
    __logBase10of2 = 3.010299956639811952137388947244930267681898814621085413104274611e-1
    if not ( type(sigfigs) is int or np.issubdtype(sigfigs, np.integer)):
        raise TypeError( "RoundToSigFigs: sigfigs must be an integer." )

    if not np.all(np.isreal( x )):
        raise TypeError( "RoundToSigFigs: all x must be real." )

    if sigfigs <= 0:
        raise ValueError( "RoundtoSigFigs: sigfigs must be positive." )

    mantissas, binaryExponents = np.frexp( x )

    decimalExponents = __logBase10of2 * binaryExponents
    intParts = np.floor(decimalExponents)

    mantissas *= 10.0**(decimalExponents - intParts)

    return np.around( mantissas, decimals=sigfigs - 1 ) * 10.0**intParts

def create_subplot(files, path, outfile):
    ''' Merge multiple plots into one image and
        index them with a letter
        
    Args:
        files: list of files to merge
        path: path to input files
        outfile: name of output file
    '''
    result = PIL.Image.new("RGB", (2000,1600), 'white')

    for index, f in enumerate(files):
        
        img = PIL.Image.open(path+f)
        
        img.thumbnail((1000, 800), PIL.Image.ANTIALIAS)
        x = index // 2 * 1000
        y = index % 2 * 800
        w, h = img.size
        print('pos {0},{1} size {2},{3}'.format(x, y, w, h))
        result.paste(img, (x, y, x + w, y + h))

    # add subfigure labels. Verdanna.tff must be in the same as this script
    fnt = PIL.ImageFont.truetype("fonts/Verdana.ttf", 30)

    PIL.ImageDraw.Draw(result).text((50,50), "a", fill=0, font=fnt)
    PIL.ImageDraw.Draw(result).text((1050,50), "b", fill=0, font=fnt)
    PIL.ImageDraw.Draw(result).text((50,780), "c", fill=0, font=fnt)
    PIL.ImageDraw.Draw(result).text((1050,780), "d", fill=0, font=fnt)

    result.save(path+outfile)

def figure_box(f, path, msg, outfile, extend=100, font="fonts/Verdana.ttf", font_size=20, x_text=0):
    ''' Add a figure box below the given image

    Args:
        f: filename of the image which will be used as input
        path: path to f
        msg: message to place in the figure box
        outfile: name of output
        extend: increase the y-axis of the image f by the given number
        font: font type to use
        font_size: font size
        x_text: specify the position where the text begins on the x-axis
    '''
    img = PIL.Image.open(path+f)
    x, y = img.size
    result = PIL.Image.new("RGB", (x, y+extend), 'white')
    result.paste(img, (0, 0))

    fnt = PIL.ImageFont.truetype(font, font_size)
    draw_word_wrap(img=result, text=msg, 
                   xpos=0+x_text, ypos=y, 
                   max_width=x-(x_text*2), font=fnt)


    result.save(path+outfile)

def draw_word_wrap(img, text, xpos=0, ypos=0, max_width=130, fill=(0, 0, 0), font=None):
    ''' Draw the given ``text`` to the x and y position of the image, using
        the minimum length word-wrapping algorithm to restrict the text to
        a pixel width of ``max_width.``

        Taken from: https://gist.github.com/atorkhov/5403562

    Args:
        img: image to draw/write upon
        text: message to draw/write
        xpos: x position to begin writing
        ypos: y position to begin writing
        max_width: maximum length on y-axis before text wrapping begins
        fill: text colour
        font: font and font size
    '''
    # font = ImageFont.truetype("font/arial.ttf", 50)
    draw = PIL.ImageDraw.Draw(img)
    text_size_x, text_size_y = draw.textsize(text, font=font)
    remaining = max_width
    space_width, space_height = draw.textsize(' ', font=font)
    # use this list as a stack, push/popping each line
    output_text = []
    # split on whitespace...    
    for word in text.split(None):
        word_width, word_height = draw.textsize(word, font=font)
        if word_width + space_width > remaining:
            output_text.append(word)
            remaining = max_width - word_width
        else:
            if not output_text:
                output_text.append(word)
            else:
                output = output_text.pop()
                output += ' %s' % word
                output_text.append(output)
            remaining = remaining - (word_width + space_width)
    for text in output_text:
        draw.text((xpos, ypos), text, font=font, fill=fill)
        ypos += text_size_y
