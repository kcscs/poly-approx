
from PySide6 import QtCore, QtGui, QtWidgets
from jsonviewer import *
import sys
# from .jsonviewer import JsonViewerWidget

SCALE_FACTOR = 1.25
RENDERED_IMAGE_PATH = "/render.png"

class PhotoViewer(QtWidgets.QGraphicsView):
    coordinatesChanged = QtCore.Signal(QtCore.QPoint)
    pixelClicked = QtCore.Signal(QtCore.QPoint)
    selected_pixel = None

    def __init__(self, parent=None):
        super().__init__(parent)
        self._zoom = 0
        self._pinned = False
        self._empty = True
        self._scene = QtWidgets.QGraphicsScene(self)
        self._photo = QtWidgets.QGraphicsPixmapItem()
        self._photo.setShapeMode(
            QtWidgets.QGraphicsPixmapItem.ShapeMode.BoundingRectShape)
        self._scene.addItem(self._photo)
        self.setScene(self._scene)
        self.setTransformationAnchor(
            QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(
            QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor(30, 30, 30)))
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)

    def hasPhoto(self):
        return not self._empty

    def resetView(self, scale=1.0):
        rect = QtCore.QRectF(self._photo.pixmap().rect())
        if not rect.isNull():
            self.setSceneRect(rect)
            if (scale := max(1, scale)) == 1:
                self._zoom = 0
            if self.hasPhoto():
                unity = self.transform().mapRect(QtCore.QRectF(0, 0, 1, 1))
                self.scale(1 / unity.width(), 1 / unity.height())
                viewrect = self.viewport().rect()
                scenerect = self.transform().mapRect(rect)
                factor = min(viewrect.width() / scenerect.width(),
                             viewrect.height() / scenerect.height()) * scale
                self.scale(factor, factor)
                if not self.zoomPinned():
                    self.centerOn(self._photo)
                self.updateCoordinates()

    def setPhoto(self, pixmap=None):
        if pixmap and not pixmap.isNull():
            self._empty = False
            self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
            self._photo.setPixmap(pixmap)
        else:
            self._empty = True
            self.setDragMode(QtWidgets.QGraphicsView.DragMode.NoDrag)
            self._photo.setPixmap(QtGui.QPixmap())
        if not (self.zoomPinned() and self.hasPhoto()):
            self._zoom = 0
        self.resetView(SCALE_FACTOR ** self._zoom)
        

    def drawForeground(self, painter, rect):
        if self.selected_pixel:
            super().drawForeground(painter, rect)
            fpos = self.selected_pixel.toPointF()
            fpos.setX(fpos.x() + 0.5)
            fpos.setY(fpos.y() + 0.5)
            painter.setPen(QtGui.QColor(0, 255, 0))  # Red border for highlight
            painter.setBrush(QtGui.QColor(0, 255, 0, 50))  # Semi-transparent red fill
            painter.drawEllipse(fpos, 2, 2)  # Draw a circle around the selected pixel


    def zoomLevel(self):
        return self._zoom

    def zoomPinned(self):
        return self._pinned

    def setZoomPinned(self, enable):
        self._pinned = bool(enable)

    def zoom(self, step):
        zoom = max(0, self._zoom + (step := int(step)))
        if zoom != self._zoom:
            self._zoom = zoom
            if self._zoom > 0:
                if step > 0:
                    factor = SCALE_FACTOR ** step
                else:
                    factor = 1 / SCALE_FACTOR ** abs(step)
                self.scale(factor, factor)
            else:
                self.resetView()

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        self.zoom(delta and delta // abs(delta))

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.resetView()

    def toggleDragMode(self):
        if self.dragMode() == QtWidgets.QGraphicsView.DragMode.ScrollHandDrag:
            self.setDragMode(QtWidgets.QGraphicsView.DragMode.NoDrag)
        elif not self._photo.pixmap().isNull():
            self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)

    def updateCoordinates(self, pos=None):
        if self._photo.isUnderMouse():
            if pos is None:
                pos = self.mapFromGlobal(QtGui.QCursor.pos())
            point = (self.mapToScene(pos) - QtCore.QPointF(0.5, 0.5)).toPoint()
        else:
            point = QtCore.QPoint()
        self.coordinatesChanged.emit(point)

    def mouseMoveEvent(self, event):
        self.updateCoordinates(event.position().toPoint())
        super().mouseMoveEvent(event)
    
    def mouseDoubleClickEvent(self, event):
        pos = event.position().toPoint()
        if self._photo.isUnderMouse():
            if pos is None:
                pos = self.mapFromGlobal(QtGui.QCursor.pos())
            point = (self.mapToScene(pos) - QtCore.QPointF(0.5, 0.5)).toPoint()
        else:
            point = QtCore.QPoint()
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            self.selected_pixel = point
            self._scene.update()
            self.pixelClicked.emit(point)
            
    def leaveEvent(self, event):
        self.coordinatesChanged.emit(QtCore.QPoint())
        super().leaveEvent(event)


class Window(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.viewer = PhotoViewer(self)
        self.viewer.coordinatesChanged.connect(self.handleCoords)
        self.viewer.pixelClicked.connect(self.selectPixel)
        self.labelCoords = QtWidgets.QLabel(self)
        self.labelCoords.setAlignment(
            QtCore.Qt.AlignmentFlag.AlignRight |
            QtCore.Qt.AlignmentFlag.AlignCenter)
        self.buttonOpen = QtWidgets.QPushButton(self)
        self.buttonOpen.setText('Open experiment folder')
        self.buttonOpen.clicked.connect(self.handleOpen)
        self.buttonPin = QtWidgets.QPushButton(self)
        self.buttonPin.setText('Pin Zoom')
        self.buttonPin.setCheckable(True)
        self.buttonPin.toggled.connect(self.viewer.setZoomPinned)
        self.jsonviewer = JsonViewerWidget()
        self.splitter = QtWidgets.QSplitter(self)
        self.splitter.addWidget(self.viewer)
        self.splitter.addWidget(self.jsonviewer)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(self.splitter, 0, 0, 1, 3)
        layout.addWidget(self.buttonOpen, 1, 0, 1, 1)
        layout.addWidget(self.buttonPin, 1, 1, 1, 1)
        layout.addWidget(self.labelCoords, 1, 2, 1, 1)
        layout.setColumnStretch(2, 2)
        self._path = "buggy_sphere"
        if len(sys.argv) > 1:
            self._path = sys.argv[1]
        self.openExperimetnFolder(self._path)

    def handleCoords(self, point):
        if not point.isNull():
            self.labelCoords.setText(f'{point.x()}, {point.y()}')
        else:
            self.labelCoords.clear()

    def handleOpen(self):
        if (start := self._path) is None:
            start = QtCore.QStandardPaths.standardLocations(
                QtCore.QStandardPaths.StandardLocation.HomeLocation)[0]
        if path := QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Open Folder', start):
            self.openExperimetnFolder(path)
            
    def openExperimetnFolder(self, path):
        self.labelCoords.clear()
        image_path = path + "/render.png"
        if not (pixmap := QtGui.QPixmap(image_path)).isNull():
            self.viewer.setPhoto(pixmap)
            self._path = path
            res = (pixmap.width(), pixmap.height())
        else:
            QtWidgets.QMessageBox.warning(self, 'Error',
                f'<br>Could not load image file:<br>'
                f'<br><b>{path}</b><br>'
                )
        metadata_path = path+"/metadata.json"
        with open(metadata_path, 'r') as file:
            json_data = json.load(file)
            self.prepareMetadata(json_data, res)

    def prepareMetadata(self, json_data, resolution):
        width, height = resolution
        self.pixelData = []
        json_pixels = json_data["pixels"]
        pixelidx = 0
        for y in range(height):
            row = []
            for x in range(width):
                row.append(json_pixels[pixelidx])
                pixelidx += 1
            self.pixelData.append(row)
        
    
    def selectPixel(self, point):
        x,y = (point.x(),point.y())
        self.jsonviewer.display_json(self.pixelData[y][x])
        print("clicked: ",point)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.setGeometry(500, 300, 800, 600)
    window.show()
    sys.exit(app.exec())
