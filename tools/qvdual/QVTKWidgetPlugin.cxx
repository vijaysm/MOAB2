/*
 * Copyright 2004 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

#include "QVTKWidgetPlugin.h"
#include "QVTKWidget.h"

#include "qpixmap.h"
#if QT_VERSION >= 0x040000
#include "qplugin.h"
#endif

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkDataSetMapper.h"
#include "vtkPolyData.h"
#include "vtkElevationFilter.h"
#include "vtkActor.h"

#include "QVTKWidget.xpm"

// macro for debug printing
#define qDebug(a)
//#define qDebug(a) printf(a)

QVTKWidgetPlugin::QVTKWidgetPlugin()
{
  qDebug("QVTKWidgetPlugin instantiated\n");
}

QVTKWidgetPlugin::~QVTKWidgetPlugin()
{
  qDebug("QVTKWidgetPlugin destructed\n");
}

QStringList QVTKWidgetPlugin::keys() const
{
  qDebug("QVTKWidgetPlugin::keys\n");
  QStringList list;
  list << "QVTKWidget";
  return list;
}

QWidget* QVTKWidgetPlugin::create( const QString& key, QWidget* parent, const char* name)
{
  qDebug("QVTKWidgetPlugin::create\n");
  if(key == "QVTKWidget")
  {
#if QT_VERSION >= 0x040000
    QVTKWidget* widget = new QVTKWidget(parent);
    widget->setObjectName(name);
#else
    QVTKWidget* widget = new QVTKWidget(parent, name);
#endif
    // gotta make a renderer so we get a nice black background in the designer
    vtkRenderer* ren = vtkRenderer::New();
    widget->GetRenderWindow()->AddRenderer(ren);

    // also for fun, let's make a cylinder and put it in the window
    // this REALLY lets the user know that a QVTKWidget works in the designer
    vtkSphereSource* cyl = vtkSphereSource::New();
    vtkElevationFilter* ele = vtkElevationFilter::New();
    ele->SetLowPoint(0.0, -0.5, 0.0);
    ele->SetHighPoint(0.0, 0.5, 0.0);
    ele->SetInput(cyl->GetOutput());
    vtkDataSetMapper* mapper = vtkDataSetMapper::New();
    mapper->SetInput(ele->GetOutput());
    ele->Delete();
    cyl->Delete();
    vtkActor* actor = vtkActor::New();
    actor->SetMapper(mapper);
    mapper->Delete();
    ren->AddProp(actor);
    actor->Delete();
    ren->Delete();

    // return the widget
    return widget;
  }
  return 0;
}

QString QVTKWidgetPlugin::group( const QString& feature ) const
{
  qDebug("QVTKWidgetPlugin::group\n");
  if(feature == "QVTKWidget")
    return "Display";
  return QString::null;
}

#if QT_VERSION < 0x040000
QIconSet QVTKWidgetPlugin::iconSet( const QString& ) const
{
  qDebug("QVTKWidgetPlugin::iconSet\n");
  return QIconSet( QPixmap( QVTKWidget_image ) );
}
#else
QIcon QVTKWidgetPlugin::iconSet( const QString& ) const
{
  qDebug("QVTKWidgetPlugin::iconSet\n");
  return QIcon( QPixmap( QVTKWidget_image ) );
}
#endif

QString QVTKWidgetPlugin::includeFile( const QString& feature ) const
{
  qDebug("QVTKWidgetPlugin::includeFile\n");
  if ( feature == "QVTKWidget" )
    return "QVTKWidget.h";
  return QString::null;
}

QString QVTKWidgetPlugin::toolTip( const QString& feature ) const
{
  qDebug("QVTKWidgetPlugin::toolTip\n");
  if(feature == "QVTKWidget")
    return "Qt VTK Widget";
  return QString::null;
}

QString QVTKWidgetPlugin::whatsThis( const QString& feature ) const
{
  qDebug("QVTKWidgetPlugin::whatsThis\n");
  if ( feature == "QVTKWidget" )
    return "A Qt/VTK Graphics Window";
  return QString::null;
}

bool QVTKWidgetPlugin::isContainer( const QString& ) const
{
  qDebug("QVTKWidgetPlugin::isContainer\n");
  return false;
}


Q_EXPORT_PLUGIN( QVTKWidgetPlugin )



