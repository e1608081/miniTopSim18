import os
import sys
import matplotlib.pyplot as plt

def plot(fname_, fname2_=None):
    """Visualizes the data given in input string fname.

        :param fname:   string containing the path to the file
    """
    #module global vars
    global fname, fname2, data1, data2, states
    fname = fname_
    fname2 = fname2_
    
    # States for line deletion and boundary
    states = {}
    states['delete'] = False
    states['boundary'] = True
    states['xup'] = 0
    states['xlo'] = 0
    states['yup'] = 0
    states['ylo'] = 0
    states['num'] = 0
    states['dL'] = 1
    # Override standard settings
    plt.rcParams['keymap.save'] = ['ctrl+s']
    plt.rcParams['keymap.all_axes'] = []
    plt.rcParams['keymap.home'] = ['h', 'home']
    plt.rcParams['keymap.quit'] = ['q']
    plt.rcParams['keymap.save'] = ['ctrl+s']
    plt.rcParams['keymap.yscale'] = []
    # Start event handler for key press
    
    # Load Surfaces
    data1 = loadFile(fname)
    data2 = None
    if fname2 is not None:
        data2 = loadFile(fname2)
    
    fig = plt.figure(1)
    cid = fig.canvas.mpl_connect('key_press_event', keyPressEventFunction)
    # Start with first surface only
    plotRewind()
    plt.show()
    # Close event handling
    fig.canvas.mpl_disconnect(cid)
    
def loadFile(fname):
    """Read topography information from file and store them for later usage.

        :param fname: string containing the path to the file
        
        :return
        retruns a list surfaces
        with every surface beeing a dictonary with
            time: simulation Time
            numel: number of elements of the Surface
            "xVals": x Values of the Surface as list
            "yVals": y Values of the Surface as list
    """
    # Load surfaces
    fHandle = open(fname, "r")
    # Place holders for surfaces, times and number of points per surface
    surfaces = []
    times = []
    numels = []
    # Create variables to store one surface
    xValues = []
    yValues = []
    for line in fHandle:
        # Check where new surface begins and get time and number of elements
        if 'surface:' in line:
            splittedLine = line.split(' ')
            time = splittedLine[1]
            numel = splittedLine[2]
            # Store extracted surface, time and number of elements
            if xValues: #do not append if list is empty; needed because of first line
                surfaces.append( (xValues, yValues) )
            
            times.append(float(time))
            numels.append(float(numel))
            # Reset lists to be empty
            xValues = []
            yValues = []
        else:
            # Extract xy pair from line
            xyValuePair = line.rstrip('\n').split(' ')
            xValue = xyValuePair[0]
            yValue = xyValuePair[1]
            xValues.append(float(xValue))
            yValues.append(float(yValue))
    # Store last line
    surfaces.append( (xValues, yValues) )
    # Close file and return extracted surfaces, times and number of points    
    fHandle.close()
    
    merged_list = []
    for i in range(len(surfaces)):
        merged_list.append({
            "time": times[i],
            "numel": numels[i],
            "xVals": surfaces[i][0],
            "yVals": surfaces[i][1],
        })
    
    return merged_list

def keyPressEventFunction(event):
    """Uses the surfaces structure and shows them to the selected settings.
        :param event:    event handle
    """
    global data1, data2, states, fname
    
    # Change boundaries
    if event.key == 'b':
        states['boundary'] = not states['boundary']
        plt.draw()

    # Step 'steps' forward
    if event.key == ' ':
        if states['num'] + 1 < len(data1):
            if states['num'] + states['dL'] < len(data1):
                states['num'] += states['dL']
            
            else: #only advance as far as possible
                states['num'] = len(data1) - 1
            
            surface1 = data1[states['num']]
            xValues1 = surface1["xVals"]
            yValues1 = surface1["yVals"]
            
            if states['delete']:
                plt.cla()
                plt.grid(True,'major')
                plt.xlabel('x-values in nm')
                plt.ylabel('y-values in nm')
            
            line, = plt.plot(xValues1, yValues1, '.-')
            
            xValues2 = None
            yValues2 = None
            if data2 is not None:
                surface2 = get_best_surface(data2, data1[states['num']]['time'])
                xValues2 = surface2["xVals"]
                yValues2 = surface2["yVals"]
                plt.plot(xValues2, yValues2, "--", color=line.get_color())
            
            # check if plot limits havo to be adopted and do so if 
            if states['boundary']:
                adoptBoundaries(xValues1, yValues1, xValues2, yValues2, states)
            
            plt.xlim(states['xlo'] - 1, states['xup'] + 1)
            plt.ylim(states['ylo'] - 1, states['yup'] + 1)
            plt.draw()


    # Change between delete state
    if event.key == 'd':
        states['delete'] = not states['delete']

    # Numbers pressed
    if str.isnumeric(event.key):
        states['dL'] = 2**int(event.key)

    # clear the figure and reset it
    if event.key == 'r':
        plotRewind()

    # change aspect ratio
    if event.key == 'a':
        if plt.axes().get_aspect() is 'equal':
            plt.axes().set_aspect( 'auto' )
        else:
            plt.axes().set_aspect( 'equal' )
        plt.draw()

    # show last line
    if event.key == 'l':
        surface = data1[-1]
        xValues = surface["xVals"]
        yValues = surface["yVals"]
        if states['delete']:
            plt.cla()
            plt.grid(True,'major')
            plt.xlabel('x-values in nm')
            plt.ylabel('y-values in nm')
        
        line, = plt.plot(xValues, yValues)
        
        xValues2 = None
        yValues2 = None
        if data2 is not None:
            surface2 = get_best_surface(data2, data1[states['num']]['time'])
            xValues2 = surface2["xVals"]
            yValues2 = surface2["yVals"]
            plt.plot(xValues2, yValues2, "--", color=line.get_color())
            
        if states['boundary']:
            adoptBoundaries(xValues, yValues, xValues2, yValues2, states)
        plt.xlim(states['xlo'] - 1, states['xup'] + 1)
        plt.ylim(states['ylo'] - 1, states['yup'] + 1)
            
        plt.draw()

    # Saving current figure
    if event.key == 's':
        plt.savefig(fname[:-4] + '.png', format='png')

def plotRewind():
    """Resets the plot to only show the first line."""
    global data1, states
    states['num'] = 0
    
    # Read values from structure
    surface1 = data1[0]
    xValues1 = surface1["xVals"]
    yValues1 = surface1["yVals"]
    
    # Plot the selected surface
    plt.cla()
    line, = plt.plot(xValues1, yValues1)
    
    xValues2 = None
    yValues2 = None
    if data2 is not None:
        surface2 = get_best_surface(data2, data1[0]['time'])
        xValues2 = surface2["xVals"]
        yValues2 = surface2["yVals"]
        plt.plot(xValues2, yValues2, "--", color=line.get_color())
        
    # check if plot limits havo to be adopted and do so if 
    if states['boundary']:
        adoptBoundaries(xValues1, yValues1, xValues2, yValues2, states)
    plt.xlim(states['xlo'] - 1, states['xup'] + 1)
    plt.ylim(states['ylo'] - 1, states['yup'] + 1)
        
    plt.grid(True,'major')
    plt.xlabel('x-values in nm')
    plt.ylabel('y-values in nm')
    plt.draw()

def adoptBoundaries(xValues, yValues, xValues2, yValues2, states):
    """Set the limits according to the given data and store them.

        :param xValues:     x coordinates of the points of the surrface
        :param yValues:     y coordinates of the points of the surrface
        :param states:      variable to store the current limits
    """
    states['xlo'] = min(xValues)
    states['xup'] = max(xValues)
    states['ylo'] = min(yValues)
    states['yup'] = max(yValues)
    
    if xValues2 is not None and yValues2 is not None:
        states['xlo'] = min(states['xlo'], min(xValues2))
        states['xup'] = max(states['xup'], max(xValues2))
        states['ylo'] = min(states['ylo'], min(yValues2))
        states['yup'] = max(states['yup'], max(yValues2))

def get_best_surface(data, time):
    """Searches for the nearest dataset for a given time"""
    best_surface = None
    smallest_diff = None
    
    for surface in data:
        diff = abs(surface['time'] - time)
        if smallest_diff is None or diff < smallest_diff:
            #we found a better one
            smallest_diff = diff
            best_surface = surface
    
    return best_surface
            

if __name__ == "__main__":
    pass
    if len(sys.argv) > 1:
        dataFile = sys.argv[1]
        if os.path.isfile(dataFile):
            plot(dataFile)
        else:
            print('No such file: ' + dataFile)
    else:
        print('No filename given. Aborted!')
