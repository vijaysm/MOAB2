/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright 2004 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/

#ifndef Q_VTK_WIDGET_H
#define Q_VTK_WIDGET_H

#include <qwidget.h>
#include <qtimer.h>

class vtkRenderWindow;
class QVTKInteractor;
#include <vtkRenderWindowInteractor.h>

#if defined(WIN32) && defined(BUILD_SHARED_LIBS)
#if defined(QVTK_EXPORTS) || defined(QVTKWidgetPlugin_EXPORTS)
#define QVTK_EXPORT __declspec( dllexport )
#else
#define QVTK_EXPORT __declspec( dllimport ) 
#endif
#else
#define QVTK_EXPORT
#endif

//! QVTKWidget displays a VTK window in a Qt window.
class QVTK_EXPORT QVTKWidget : public QWidget
{
  Q_OBJECT

  public:
#if QT_VERSION < 0x040000
    //! constructor for Qt 3
    QVTKWidget(QWidget* parent = NULL, const char* name = NULL, Qt::WFlags f = 0);
#else
    //! constructor for Qt 4
    QVTKWidget(QWidget* parent = NULL, Qt::WFlags f = 0);
#endif
    //! destructor
    virtual ~QVTKWidget();

    //! set the vtk render window, if you wish to use your own vtkRenderWindow
    void SetRenderWindow(vtkRenderWindow*);
    
    //! get the vtk render window
    vtkRenderWindow* GetRenderWindow();
    
    //! get the Qt/vtk interactor that was either created by default or set by the user
    QVTKInteractor* GetInteractor();
    
    //! handle reparenting of this widget
#if QT_VERSION < 0x040000
    virtual void reparent(QWidget* parent, Qt::WFlags f, const QPoint& p, bool showit);
#else
    virtual void setParent(QWidget* parent, Qt::WFlags f);
#endif
    
    //! handle hides
    virtual void hide();
    //! handle shows
    virtual void show();
  
  protected:
    //! overloaded resize handler
    virtual void resizeEvent(QResizeEvent* event);
    //! overloaded move handler
    virtual void moveEvent(QMoveEvent* event);
    //! overloaded paint handler
    virtual void paintEvent(QPaintEvent* event);

    //! overloaded mouse press handler
    virtual void mousePressEvent(QMouseEvent* event);
    //! overloaded mouse move handler
    virtual void mouseMoveEvent(QMouseEvent* event);
    //! overloaded mouse release handler
    virtual void mouseReleaseEvent(QMouseEvent* event);
    //! overloaded key press handler
    virtual void keyPressEvent(QKeyEvent* event);
    //! overloaded key release handler
    virtual void keyReleaseEvent(QKeyEvent* event);
    //! overloaded enter event
    virtual void enterEvent(QEvent*);
    //! overloaded leave event
    virtual void leaveEvent(QEvent*);
#ifndef QT_NO_WHEELEVENT
    //! overload wheel mouse event
    virtual void wheelEvent(QWheelEvent*);
#endif
    //! overload focus event
    virtual void focusInEvent(QFocusEvent*);
    //! overload focus event
    virtual void focusOutEvent(QFocusEvent*);
    //! overload Qt's event() to capture more keys
    bool event( QEvent* e );
    
    //! the vtk render window
    vtkRenderWindow* mRenWin;

#if defined(Q_WS_MAC)
    void macFixRect();
    virtual void setRegionDirty(bool);
    virtual void macWidgetChangedWindow();
#endif
  private slots:
    void internalMacFixRect();

  private:
    //! unimplemented operator=
    QVTKWidget const& operator=(QVTKWidget const&);
    //! unimplemented copy
    QVTKWidget(const QVTKWidget&);

};


//! Qt/VTK interactor class
class QVTK_EXPORT QVTKInteractor : public QObject, public vtkRenderWindowInteractor
{
  Q_OBJECT
public:
  //! allocation method
  static QVTKInteractor* New();
  vtkTypeMacro(QVTKInteractor,vtkRenderWindowInteractor);

  //! overloaded terminiate app
  virtual void TerminateApp();
  //! overloaded start method
  virtual void Start();
  //! overloaded create timer method
  virtual int CreateTimer(int);
  //! overloaded destroy timer method
  virtual int DestroyTimer();

public slots:
  //! timer event slot
  virtual void TimerEvent();

protected:
  //! constructor
  QVTKInteractor();
  //! destructor
  ~QVTKInteractor();
private:

  //! timer
  QTimer mTimer;
  
  //! unimplemented copy
  QVTKInteractor(const QVTKInteractor&);
  //! unimplemented operator=
  void operator=(const QVTKInteractor&);

};


#endif


