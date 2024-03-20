import QtQuick 2.12
import QtQuick.Window 2.12

Window {
    id: window
    width: 640
    height: 480
    visible: true
    title: qsTr("Hello World")

    GridMenu {
        anchors.fill: parent

        menuItems: [
            {name: "Option 1", icon: "icon1.png", targetItem: "Page1.qml"},
            {name: "Option 2", icon: "icon2.png", targetItem: "Page2.qml"},
            {name: "Option 3", icon: "icon3.png", targetItem: "Page3.qml"},
            {name: "Option 4", icon: "icon4.png", targetItem: "Page4.qml"},
            {name: "Option 1", icon: "icon1.png", targetItem: "Page1.qml"},
            {name: "Option 2", icon: "icon2.png", targetItem: "Page2.qml"},
            {name: "Option 3", icon: "icon3.png", targetItem: "Page3.qml"},
            {name: "Option 4", icon: "icon4.png", targetItem: "Page4.qml"},
            {name: "Option 1", icon: "icon1.png", targetItem: "Page1.qml"},
            {name: "Option 2", icon: "icon2.png", targetItem: "Page2.qml"},
            {name: "Option 3", icon: "icon3.png", targetItem: "Page3.qml"},
            {name: "Option 4", icon: "icon4.png", targetItem: "Page4.qml"},
            {name: "Option 1", icon: "icon1.png", targetItem: "Page1.qml"},
            {name: "Option 2", icon: "icon2.png", targetItem: "Page2.qml"},
            {name: "Option 3", icon: "icon3.png", targetItem: "Page3.qml"},
            {name: "Option 4", icon: "icon4.png", targetItem: "Page4.qml"},
        ]
    }
}
