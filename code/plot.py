import os
import sys
import copy
import matplotlib.pyplot as plt

def plot(fname):
    """Visualizes the data given in input string fname.

        :param fname:   string containing the path to the file
    """
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
    # Load Surfaces
    times, numels, surfaces = loadFile(fname)
    # Override standard settings
    plt.rcParams['keymap.save'] = ['ctrl+s']
    plt.rcParams['keymap.all_axes'] = []
    plt.rcParams['keymap.home'] = ['h', 'home']
    plt.rcParams['keymap.quit'] = ['q']
    plt.rcParams['keymap.save'] = ['ctrl+s']
    plt.rcParams['keymap.yscale'] = []
    # Start event handler for key press
    fig = plt.figure(1)
    cid = fig.canvas.mpl_connect('key_press_event', lambda event: keyPressEventFunction(event, times, numels, surfaces, states, fname))
    # Start with first surface only
    plotRewind(times, numels, surfaces, states)
    plt.show()
    # Close event handling
    fig.canvas.mpl_disconnect(cid)
    
def loadFile(fname):
    """Read topography information from file and store them for later usage.

        :param fname: string containing the path to the file
    """
    # Load surfaces
    fHandle = open( fname, "r" )
    # Place holders for surfaces, times and number of points per surface
    surfaces = list()
    times = list()
    numels = list()
    # Create variables to store one surface
    xValues = list()
    yValues = list()
    for line in fHandle:
        # Check where new surface begins and get time and number of elements
        if 'surface:' in line:
            splittedLine = line.split(' ')
            time = splittedLine[1]
            numel = splittedLine[2]
            # Store extracted surface, time and number of elements
            if xValues:
                xyValuePairSorted = zip(*zip(xValues, yValues))
                surfaces.append( xyValuePairSorted )
            
            times.append(float(time))
            numels.append(float(numel))
            # Reset lists to be empty
            xValues.clear()
            yValues.clear()
        else:
            # Extract xy pair from line
            xyValuePair = line.rstrip('\n').split(' ')
            xValue = xyValuePair[0]
            yValue = xyValuePair[1]
            xValues.append(float(xValue))
            yValues.append(float(yValue))
    # Store last line
    xyValuePairSorted = zip(*zip(xValues, yValues))
    surfaces.append(xyValuePairSorted)
    # Close file and return extracted surfaces, times and number of points    
    fHandle.close()
    return times, numels, surfaces

def keyPressEventFunction(event, times, numels, surfaces, states, fname):
    """Uses the surfaces structure and shows them to the selected settings.

        :param event:    event handle
        :param times:    time vector extracted from file
        :param numels:   number of elements in a surface
        :param surfaces: list of zipped lists holding the position
        :param states:   settings for plots, delete last lines and adopt
                         limits of axis, current axis limits
        :param fname:    string file name
    """
    surfacesCopy = copy.deepcopy( surfaces )
    # Change boundaries
    if event.key == 'b':
        states['boundary'] = not states['boundary']
        plt.draw()

    # Step 'steps' forward
    if event.key == ' ':
        if states['num'] + states['dL'] < len(surfacesCopy):
            states['num'] += states['dL']
            surface = copy.copy( surfacesCopy[states['num']] )
            xValues, yValues = surface
            
            if states['delete']:
                plt.cla()
            if states['boundary']:
                adoptBoundaries(xValues, yValues, states)

            plt.plot(xValues, yValues)
            plt.xlim( states['xlo'] - 1, states['xup'] + 1 )
            plt.ylim( states['ylo'] - 1, states['yup'] + 1 )
            plt.draw()

            


    # Change between delete state
    if event.key == 'd':
        states['delete'] = not states['delete']

    # Numbers pressed
    if str.isnumeric(event.key):
        states['dL'] = 2**int(event.key)

    # clear the figure and reset it
    if event.key == 'r':
        plotRewind(times, numels, surfacesCopy, states)

    # change aspect ratio
    if event.key == 'a':
        if plt.axes().get_aspect() is 'equal':
            plt.axes().set_aspect( 'auto' )
        else:
            plt.axes().set_aspect( 'equal' )
        plt.draw()

    # show last line
    if event.key == 'l':
        surface = surfacesCopy[-1]
        xValues, yValues = surface
        if states['delete']:
            plt.cla()
        if states['boundary']:
            adoptBoundaries(xValues, yValues, states)
        plt.xlim( states['xlo'] - 1, states['xup'] + 1 )
        plt.ylim( states['ylo'] - 1, states['yup'] + 1 )
        plt.plot(xValues, yValues)
        plt.draw()

    # Saving current figure
    if event.key == 's':
        plt.savefig( fname[:-4] + '.png', format='png' )

def plotRewind(times, numels, surfaces, states):
    """Resets the plot to only show the first line.

        :param event:   event handler variable
        :param times:   time variable for each line
        :param numels:  number of points for the lines
        :param surface: coordinate values of the surface
        :param states:  variable storing the delete and boundary state for
                        the plot and the current x and y limits of the plot
    """
    states['num'] = 0
    # Read values from structure
    surfacesCopy = copy.deepcopy( surfaces )
    surface = surfacesCopy[1]
    xValues, yValues = surface
    # check if plot limits havo to be adopted and do so if 
    if states['boundary']:
        adoptBoundaries(xValues, yValues, states)
    # Plot the selected surface
    plt.cla()
    plt.xlim( states['xlo'] - 1, states['xup'] + 1 )
    plt.ylim( states['ylo'] - 1, states['yup'] + 1 )
    plt.plot(xValues, yValues)
    plt.grid(True,'major')
    plt.xlabel('x-values in nm')
    plt.ylabel('y-values in nm')
    plt.draw()

def adoptBoundaries(xValues, yValues, states):
    """Set the limits according to the given data and store them.

        :param xValues:     x coordinates of the points of the surrface
        :param yValues:     y coordinates of the points of the surrface
        :param states:      variable to store the current limits
    """
    states['xlo'] = xValues[0]
    states['xup'] = xValues[-1]
    states['ylo'] = min(yValues)
    states['yup'] = max(yValues)






if __name__ == "__main__":
    pass
    if len(sys.argv) > 1:
        dataFile = sys.argv[1]
        if os.path.isfile(dataFile):
            plot( dataFile )
        else:
            print('No such file: ' + dataFile)
    else:
        print('No filename given. Aborted!')
