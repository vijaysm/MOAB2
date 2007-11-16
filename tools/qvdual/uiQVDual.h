/****************************************************************************
** Form interface generated from reading ui file 'uiQVDual.ui'
**
** Created: Mon Nov 5 08:55:16 2007
**      by: The User Interface Compiler ($Id: qt/main.cpp   3.3.7   edited Aug 31 2005 $)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#ifndef UIQVDUAL_H
#define UIQVDUAL_H

#include <qvariant.h>
#include <qmainwindow.h>
#include <map>
#include <set>
#include "MBInterface.hpp"

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class QAction;
class QActionGroup;
class QToolBar;
class QPopupMenu;
class QListView;
class QListViewItem;
class QPushButton;
class QLabel;
class QLineEdit;
class QVTKWidget;
class vtkUnstructuredGrid;
class vtkRenderer;
class vtkRenderWindowInteractor;
class vtkProperty;
class CropToolPopup;
class DrawDual;

class uiQVDual : public QMainWindow
{
    Q_OBJECT

public:
    uiQVDual( QWidget* parent = 0, const char* name = 0, WFlags fl = WType_TopLevel );
    ~uiQVDual();

    QListView* ActorListView1;
    QListView* TagListView1;
    QPushButton* APbutton;
    QPushButton* FOCbutton;
    QPushButton* FSbutton;
    QPushButton* Create_Dual_Button;
    QPushButton* DebugButton;
    QPushButton* negAPbutton;
    QPushButton* negFOCbutton;
    QPushButton* negFCbutton;
    QPushButton* CropTool;
    QPushButton* CropToolButton;
    QLabel* picklabel1;
    QLineEdit* pickline1;
    QLabel* picklabel2;
    QLineEdit* pickline2;
    QVTKWidget* vtkWidget;
    QMenuBar *MenuBar;
    QPopupMenu *fileMenu;
    QPopupMenu *helpMenu;
    QPopupMenu *Window;
    QPopupMenu *Display;
    QAction* fileNewAction;
    QAction* fileOpenAction;
    QAction* fileSaveAction;
    QAction* fileSaveAsAction;
    QAction* filePrintAction;
    QAction* fileExitAction;
    QAction* helpContentsAction;
    QAction* helpIndexAction;
    QAction* helpAboutAction;
    QAction* windowTagViewAction;
    QAction* windowSheetViewAction;
    QAction* windowChordViewAction;
    QAction* windowActorViewAction;
    QAction* displayVisibleAction;
    QAction* displayDrawAction;
    QAction* displayWireframeShadedAction;
    QAction* displayInvertSelectionAction;
    QAction* displayInvisibleAction;
    QAction* displayDisplayAllActions;
    QAction* displayDrawSheetAction;

    bool computeDual;

public slots:
    virtual void fileNew();
    virtual void fileOpen();
    virtual void fileOpen( const QString & filename );
    virtual void fileSave();
    virtual void fileSaveAs( const QString & filename );
    virtual void filePrint();
    virtual void fileExit();
    virtual void helpIndex();
    virtual void helpContents();
    virtual void helpAbout();
    virtual void constructDual();
    virtual void updateMesh();
    virtual void DebugButton_pressed();
    virtual void ActorTreeView_selectionChanged();
    virtual void TagTreeView_selectionChanged();
    virtual void displayVisible();
    virtual void displayDraw();
    virtual void displayWireframeShaded();
    virtual void displayInvertSelection();
    virtual void ActorListView1_rightButtonPressed( QListViewItem * item, const QPoint &, int );
    virtual void displayInvisible();
    virtual void TagListView1_rightButtonPressed( QListViewItem *, const QPoint &, int );
    virtual void displayDisplayAll();
    virtual void CropToolButton_clicked();
    virtual void displayDrawSheetAction_activated();
    virtual void APbutton_clicked();
    virtual void negAPbutton_clicked();
    virtual void FOCbutton_clicked();
    virtual void FSbutton_clicked();
    virtual void negFCbutton_clicked();
    virtual void fileSaveAs();
    virtual void resetDisplay();
    virtual void redrawDisplay();
    virtual void pickline1_returnPressed();

signals:
    void toggled();

protected:
    CropToolPopup *cropToolPopup;
    std::set<QListViewItem*> itemSelList;
    DrawDual *drawDual;

    QVBoxLayout* layout13;
    QHBoxLayout* layout6;
    QHBoxLayout* layout7;

protected slots:
    virtual void languageChange();

private:
    QString lastOpened;
    std::map<QListViewItem*,MBEntityHandle> itemSetMap;
    vtkUnstructuredGrid *myUG;
    int currentWin;

    void init();
    void destroy();
    virtual void updateTagList();
    virtual void updateActorList();
    virtual void updateActorContainsList( QListViewItem * item, MBEntityHandle set_handle );
    virtual void updateActorParentList( QListViewItem * item, MBEntityHandle set_handle );
    virtual void changeSetProperty( std::set<QListViewItem *> & high_sets, std::set<QListViewItem *> & unhigh_sets );
    virtual void evalItem( QListViewItem * item, const bool high, MBRange & high_mbsets, MBRange & unhigh_mbsets );
    virtual void getSelected( QListView * listv, std::set<QListViewItem *> & selected, std::set<QListViewItem *> & unselected );
    virtual void getSelected( QListView * listv, MBRange & selected, MBRange & unselected );
    virtual void getItemSets( std::set<QListViewItem *> & items, MBRange & sets );

};

#endif // UIQVDUAL_H
