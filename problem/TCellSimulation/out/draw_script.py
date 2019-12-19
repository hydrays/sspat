#### import the simple module from the paraview
import sys
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#baseDIR = '/home/hydra/work/Research/Tcell/sspat/problem/TCellSimulation/out/'
baseDIR = sys.argv[1]

# create a new 'CSV Reader'
chemokinecsv = CSVReader(FileName=[baseDIR + 'chemokine.csv'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [930, 548]

# get layout
viewLayout1 = GetLayout()

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=chemokinecsv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'
tableToPoints1.KeepAllDataArrays = 1

# create a new 'Delaunay 2D'
delaunay2D1 = Delaunay2D(Input=tableToPoints1)

# set active view
#SetActiveView(renderView1)

# set active source
#SetActiveSource(delaunay2D1)

# show data in view
delaunay2D1Display = Show(delaunay2D1, renderView1)
# trace defaults for the display properties.
delaunay2D1Display.ColorArrayName = [None, '']
delaunay2D1Display.GlyphType = 'Arrow'

# reset view to fit data
renderView1.ResetCamera()

# change solid color
# normal signal
#delaunay2D1Display.DiffuseColor = [0.6862745098039216, 0.8627450980392157, 0.9098039215686274]
#delaunay2D1Display.DiffuseColor = [0.7843137254901961, 0.9215686274509803, 0.9607843137254902]
#delaunay2D1Display.DiffuseColor = [0.803921568627451, 0.9411764705882353, 0.9803921568627451]

# weak signal
#delaunay2D1Display.DiffuseColor = [0.9215686274509803, 0.9607843137254902, 1.0]
#delaunay2D1Display.DiffuseColor = [0.9372549019607843, 0.9764705882352941, 1.0]
 
# no signal
#delaunay2D1Display.DiffuseColor = [0.9607843137254902, 0.9607843137254902, 0.9607843137254902]
delaunay2D1Display.DiffuseColor = [0.9803921568627451, 0.9803921568627451, 0.9803921568627451]

# get color transfer function/color map for 'z'
#zLUT = GetColorTransferFunction('z')

# get opacity transfer function/opacity map for 'z'
#zPWF = GetOpacityTransferFunction('z')

# create a new 'CSV Reader'
fBcellcsv = CSVReader(FileName=[baseDIR + 'FBcell.csv'])

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(Input=fBcellcsv)
tableToPoints2.XColumn = 'x'
tableToPoints2.YColumn = 'y'
tableToPoints2.ZColumn = 'z'
tableToPoints2.KeepAllDataArrays = 1

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints2)

# show data in view
tableToPoints2Display = Show(tableToPoints2, renderView1)
# trace defaults for the display properties.
tableToPoints2Display.ColorArrayName = [None, '']
tableToPoints2Display.GlyphType = 'Arrow'

# set scalar coloring
ColorBy(tableToPoints2Display, ('POINTS', 't'))

renderView1.ResetCamera()

# rescale color and/or opacity maps used to include current data range
tableToPoints2Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
#tableToPoints2Display.SetScalarBarVisibility(renderView1, False)

# get color transfer function/color map for 't'
tLUT = GetColorTransferFunction('t')

# get opacity transfer function/opacity map for 't'
#tPWF = GetOpacityTransferFunction('t')

# Properties modified on tLUT
#tLUT.ColorSpace = 'RGB'

# Properties modified on tLUT
#tLUT.RGBPoints = [1.0, 0.3333333333333333, 0.0, 0.0, 2.0, 0.12156862745098039, 0.47058823529411764, 0.7137254901960784]
#tLUT.RGBPoints = [1.0, 1, 1, 0.2, 2.0, 0.9725490196078431, 0.4627450980392157, 0.42745098039215684]

## show data in view
#tableToPoints2Display = Show(tableToPoints2, renderView1)
## trace defaults for the display properties.
#tableToPoints2Display.ColorArrayName = [None, '']
#tableToPoints2Display.GlyphType = 'Arrow'

## change solid color
#tableToPoints2Display.DiffuseColor = [0.12156862745098039, 0.47058823529411764, 0.7137254901960784]

# create a new 'CSV Reader'
domaincsv = CSVReader(FileName=[baseDIR + 'domain.csv'])

# create a new 'Table To Points'
tableToPoints3 = TableToPoints(Input=domaincsv)
tableToPoints3.XColumn = 'x'
tableToPoints3.YColumn = 'y'
tableToPoints3.ZColumn = 'z'

# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(Input=tableToPoints3)

# show data in view
delaunay3D1Display = Show(delaunay3D1, renderView1)
# trace defaults for the display properties.
delaunay3D1Display.ColorArrayName = [None, '']
delaunay3D1Display.GlyphType = 'Arrow'
delaunay3D1Display.ScalarOpacityUnitDistance = 201.12866189986227

# Properties modified on delaunay3D1Display
delaunay3D1Display.Opacity = 0.05

# create a new 'CSV Reader'
tMcellcsv = CSVReader(FileName=[baseDIR + 'TMcell.csv'])

# create a new 'Table To Points'
tableToPoints4 = TableToPoints(Input=tMcellcsv)
tableToPoints4.XColumn = 'x'
tableToPoints4.YColumn = 'y'
tableToPoints4.ZColumn = 'z'
tableToPoints4.KeepAllDataArrays = 1

# show data in view
tableToPoints4Display = Show(tableToPoints4, renderView1)
# trace defaults for the display properties.
tableToPoints4Display.ColorArrayName = [None, '']
tableToPoints4Display.GlyphType = 'Arrow'

# set scalar coloring
ColorBy(tableToPoints4Display, ('POINTS', 't'))

# rescale color and/or opacity maps used to include current data range
tableToPoints4Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
tableToPoints4Display.SetScalarBarVisibility(renderView1, False)

# get color transfer function/color map for 't'
tLUT = GetColorTransferFunction('t')

# get opacity transfer function/opacity map for 't'
# tPWF = GetOpacityTransferFunction('t')

# Properties modified on tLUT
tLUT.RGBPoints = [1.0, 0.8, 0.8, 0.8, 2.0, 0.865003, 0.865003, 0.865003, 3.0, 0.9725490196078431, 0.4627450980392157, 0.42745098039215684]

tableToPoints4Display.Opacity = 0.7

# create a new 'CSV Reader'
inactiveTcellcsv = CSVReader(FileName=[baseDIR + 'inactiveTcell.csv'])

# create a new 'Table To Points'
tableToPoints5 = TableToPoints(Input=inactiveTcellcsv)
tableToPoints5.XColumn = 'x'
tableToPoints5.YColumn = 'y'
tableToPoints5.ZColumn = 'z'

# show data in view
tableToPoints5Display = Show(tableToPoints5, renderView1)
# trace defaults for the display properties.
tableToPoints5Display.ColorArrayName = [None, '']
tableToPoints5Display.GlyphType = 'Arrow'

# Properties modified on tableToPoints1Display
tableToPoints5Display.Opacity = 0.7

# change solid color
#tableToPoints1Display.DiffuseColor = [0.705882, 0.0156863, 0.14902]# change solid color
tableToPoints5Display.DiffuseColor = [0.9725490196078431, 0.4627450980392157, 0.42745098039215684]

#### saving camera placements for all active views

# current camera placement for renderView1
#renderView1.CameraPosition = [-500, -190, 368]
#renderView1.CameraFocalPoint = [256, 256, -1.67]
#renderView1.CameraViewUp = [0.35, 0.18, 0.92]
#renderView1.CameraParallelScale = 361.410368898389

# current camera placement for renderView1
renderView1.CameraPosition = [-528.076345693851, -153.193616227026, 358.899186383466]
renderView1.CameraFocalPoint = [228.84446987953055, 294.2409529138935, -10.571272252883812]
renderView1.CameraViewUp = [0.34017899229008414, 0.18573151220164705, 0.9218362428217891]
renderView1.CameraParallelScale = 361.410368898389

renderView1.Background = [1.0, 1.0, 1.0]

renderView1.OrientationAxesVisibility = 0
# save screenshot
SaveScreenshot(baseDIR + 'pvoutput.png', magnification=1, quality=100, view=renderView1)

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
