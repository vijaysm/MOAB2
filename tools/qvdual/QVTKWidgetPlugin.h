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

#ifndef Q_VTK_WIDGET_PLUGIN_H
#define Q_VTK_WIDGET_PLUGIN_H

/****************

Plugin code to allow QVTKWidget plugin in the Qt designer

Build this and put it in your plugin path.
That could be in $QTDIR/lib/plugin or any path/plugin pointed to by qtconfig.

*************/

#include <qwidgetplugin.h>

// derive from QWidgetPlugin and implement the plugin interface
class QVTKWidgetPlugin : public QWidgetPlugin
{
  public:
    QVTKWidgetPlugin();
    ~QVTKWidgetPlugin();
    
    //! return a list of keys for what widgets this plugin makes
    QStringList keys() const;
    //! create a widget by key
    QWidget* create( const QString& key, QWidget* parent = 0, const char* name = 0);
    //! what group this plugin shows up in the designer
    QString group( const QString& ) const;
    //! the icons for the widgets
#if QT_VERSION < 0x040000
    QIconSet iconSet( const QString& ) const;
#else
    QIcon iconSet( const QString& ) const;
#endif
    //! the name of the include file for building an app with a widget
    QString includeFile( const QString& ) const;
    //! tool tip text
    QString toolTip( const QString& ) const;
    //! what's this text
    QString whatsThis( const QString& ) const;
    //! returns whether widget is a container
    bool isContainer( const QString& ) const;
};

#endif


