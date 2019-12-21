#### import the simple module from the paraview
import sys
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

baseDIR = sys.argv[1]
##baseDIR = './'

# create a new 'CSV Reader'
chemokinecsv = CSVReader(FileName=[baseDIR + 'chemokine.csv'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [930, 548]

# get layout
viewLayout1 = GetLayout()

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# place view in the layout
viewLayout1.AssignView(2, spreadSheetView1)

# show data in view
chemokinecsvDisplay = Show(chemokinecsv, spreadSheetView1)
# trace defaults for the display properties.
chemokinecsvDisplay.FieldAssociation = 'Row Data'

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=chemokinecsv)
tableToPoints1.XColumn = 'h'
tableToPoints1.YColumn = 'h'
tableToPoints1.ZColumn = 'h'

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'
tableToPoints1.KeepAllDataArrays = 1

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# hide data in view
Hide(chemokinecsv, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
viewLayout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# create a new 'Delaunay 2D'
delaunay2D1 = Delaunay2D(Input=tableToPoints1)

# show data in view
delaunay2D1Display = Show(delaunay2D1, renderView1)
# trace defaults for the display properties.
delaunay2D1Display.ColorArrayName = [None, '']
delaunay2D1Display.GlyphType = 'Arrow'

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(delaunay2D1Display, ('POINTS', 'floor'))

# rescale color and/or opacity maps used to include current data range
delaunay2D1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
delaunay2D1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 't'
tLUT = GetColorTransferFunction('floor')

# get opacity transfer function/opacity map for 't'
tPWF = GetOpacityTransferFunction('floor')

# hide color bar/color legend
delaunay2D1Display.SetScalarBarVisibility(renderView1, False)

#tLUT.RGBPoints = [1.0, 0.8, 0.8, 0.8, 2.0, 0.865003, 0.865003, 0.865003, 3.0, 0.9725490196078431, 0.4627450980392157, 0.42745098039215684]
tLUT.RGBPoints = [1.0, 0.8627450980392157, 1.0, 0.9411764705882353, 2.0, 0.7058823529411765, 0.8823529411764706, 1.0]


# Final version
# strong signal
#delaunay2D1Display.DiffuseColor = [0.7058823529411765, 0.8823529411764706, 1.0]
# weak signal
##delaunay2D1Display.DiffuseColor = [0.8627450980392157, 1.0, 0.9411764705882353]
# no signal
#delaunay2D1Display.DiffuseColor = [0.9803921568627451, 0.9803921568627451, 0.9803921568627451]

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
#delaunay2D1Display.DiffuseColor = [0.9803921568627451, 0.9803921568627451, 0.9803921568627451]

# create a new 'CSV Reader'
tMcellcsv = CSVReader(FileName=[baseDIR + 'TMcell.csv'])

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# place view in the layout
viewLayout1.AssignView(2, spreadSheetView1)

# show data in view
tMcellcsvDisplay = Show(tMcellcsv, spreadSheetView1)
# trace defaults for the display properties.
tMcellcsvDisplay.FieldAssociation = 'Row Data'

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(Input=tMcellcsv)
tableToPoints2.XColumn = 't'
tableToPoints2.YColumn = 't'
tableToPoints2.ZColumn = 't'

# Properties modified on tableToPoints2
tableToPoints2.XColumn = 'x'
tableToPoints2.YColumn = 'y'
tableToPoints2.ZColumn = 'z'
tableToPoints2.KeepAllDataArrays = 1

# show data in view
tableToPoints2Display = Show(tableToPoints2, spreadSheetView1)

# hide data in view
Hide(tMcellcsv, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
viewLayout1.Collapse(2)

renderView1.ResetCamera()

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints2)

# show data in view
tableToPoints2Display = Show(tableToPoints2, renderView1)
# trace defaults for the display properties.
tableToPoints2Display.ColorArrayName = [None, '']
tableToPoints2Display.GlyphType = 'Arrow'

# Properties modified on tableToPoints2Display
tableToPoints2Display.Opacity = 0.7

# set scalar coloring
ColorBy(tableToPoints2Display, ('POINTS', 't'))

# rescale color and/or opacity maps used to include current data range
tableToPoints2Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
tableToPoints2Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 't'
tLUT = GetColorTransferFunction('t')

# get opacity transfer function/opacity map for 't'
tPWF = GetOpacityTransferFunction('t')

# hide color bar/color legend
tableToPoints2Display.SetScalarBarVisibility(renderView1, False)

#tLUT.RGBPoints = [1.0, 0.8, 0.8, 0.8, 2.0, 0.865003, 0.865003, 0.865003, 3.0, 0.9725490196078431, 0.4627450980392157, 0.42745098039215684]
tLUT.RGBPoints = [1.0, 0.8, 0.8, 0.8, 2.0, 0.865003, 0.865003, 0.865003, 3.0, 0.9725490196078431, 0.4627450980392157, 0.42745098039215684]

tableToPoints2Display.Opacity = 0.7

# create a new 'CSV Reader'
fBcellcsv = CSVReader(FileName=[baseDIR + 'FBcell.csv'])

print "haha"
print fBcellcsv

#if fBcellcsv.NumofRows > 0
# create a new 'Table To Points'
tableToPoints3 = TableToPoints(Input=fBcellcsv)
tableToPoints3.XColumn = 'x'
tableToPoints3.YColumn = 'y'
tableToPoints3.ZColumn = 'z'
tableToPoints3.KeepAllDataArrays = 1

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints3)

# show data in view
tableToPoints3Display = Show(tableToPoints3, renderView1)
# trace defaults for the display properties.
tableToPoints3Display.ColorArrayName = [None, '']
tableToPoints3Display.GlyphType = 'Arrow'

# set scalar coloring
ColorBy(tableToPoints3Display, ('POINTS', 'h'))

renderView1.ResetCamera()

# rescale color and/or opacity maps used to include current data range
tableToPoints3Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
tableToPoints2Display.SetScalarBarVisibility(renderView1, False)

# get color transfer function/color map for 't'
tLUT = GetColorTransferFunction('h')

# get opacity transfer function/opacity map for 't'
tPWF = GetOpacityTransferFunction('h')

# Properties modified on tLUT
#tLUT.ColorSpace = 'RGB'

# Properties modified on tLUT
tLUT.RGBPoints = [1.0, 0.12156862745098039, 0.7137254901960784, 0.47058823529411764, 2.0, 0.12156862745098039, 0.47058823529411764, 0.7137254901960784]
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

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# place view in the layout
viewLayout1.AssignView(2, spreadSheetView1)

# show data in view
domaincsvDisplay = Show(domaincsv, spreadSheetView1)
# trace defaults for the display properties.
domaincsvDisplay.FieldAssociation = 'Row Data'

# create a new 'Table To Points'
tableToPoints4 = TableToPoints(Input=domaincsv)
tableToPoints4.XColumn = 'x'
tableToPoints4.YColumn = 'x'
tableToPoints4.ZColumn = 'x'

# Properties modified on tableToPoints4
tableToPoints4.YColumn = 'y'
tableToPoints4.ZColumn = 'z'

# show data in view
tableToPoints4Display = Show(tableToPoints4, spreadSheetView1)

# hide data in view
Hide(domaincsv, spreadSheetView1)

# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(Input=tableToPoints4)

# show data in view
delaunay3D1Display = Show(delaunay3D1, spreadSheetView1)

# hide data in view
Hide(tableToPoints4, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
viewLayout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(delaunay3D1)

# show data in view
delaunay3D1Display = Show(delaunay3D1, renderView1)
# trace defaults for the display properties.
delaunay3D1Display.ColorArrayName = [None, '']
delaunay3D1Display.GlyphType = 'Arrow'
delaunay3D1Display.ScalarOpacityUnitDistance = 304.2360599886002

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

delaunay3D1Display.Opacity = 0.05

# create a new 'CSV Reader'
inactiveTcellcsv = CSVReader(FileName=[baseDIR + 'inactiveTcell.csv'])

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# place view in the layout
viewLayout1.AssignView(2, spreadSheetView1)

# show data in view
inactiveTcellcsvDisplay = Show(inactiveTcellcsv, spreadSheetView1)
# trace defaults for the display properties.
inactiveTcellcsvDisplay.FieldAssociation = 'Row Data'

# create a new 'Table To Points'
tableToPoints5 = TableToPoints(Input=inactiveTcellcsv)
tableToPoints5.XColumn = 'x'
tableToPoints5.YColumn = 'x'
tableToPoints5.ZColumn = 'x'

# Properties modified on tableToPoints5
tableToPoints5.YColumn = 'y'
tableToPoints5.ZColumn = 'z'

# show data in view
tableToPoints5Display = Show(tableToPoints5, spreadSheetView1)

# hide data in view
Hide(inactiveTcellcsv, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
viewLayout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints5)

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
renderView1.CameraPosition = [1190.4042920770053, -640.6984507019484, 505.685075109063]
renderView1.CameraFocalPoint = [513.3297866233299, 539.482098990153, 4.750712460976326]
renderView1.CameraViewUp = [-0.14552485054532527, 0.314645909724804, 0.9379874569349056]
renderView1.CameraParallelScale = 549.4182526985146

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

renderView1.OrientationAxesVisibility = 0
# save screenshot
SaveScreenshot(baseDIR + 'pvoutput.png', magnification=1, quality=100, view=renderView1)
