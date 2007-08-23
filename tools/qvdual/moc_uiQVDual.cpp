/****************************************************************************
** uiQVDual meta object code from reading C++ file 'uiQVDual.h'
**
** Created: Wed Aug 22 09:12:40 2007
**      by: The Qt MOC ($Id: qt/moc_yacc.cpp   3.3.7   edited Oct 19 16:22 $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "uiQVDual.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *uiQVDual::className() const
{
    return "uiQVDual";
}

QMetaObject *uiQVDual::metaObj = 0;
static QMetaObjectCleanUp cleanUp_uiQVDual( "uiQVDual", &uiQVDual::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString uiQVDual::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "uiQVDual", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString uiQVDual::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "uiQVDual", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* uiQVDual::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QMainWindow::staticMetaObject();
    static const QUMethod slot_0 = {"fileNew", 0, 0 };
    static const QUMethod slot_1 = {"fileOpen", 0, 0 };
    static const QUParameter param_slot_2[] = {
	{ "filename", &static_QUType_QString, 0, QUParameter::In }
    };
    static const QUMethod slot_2 = {"fileOpen", 1, param_slot_2 };
    static const QUMethod slot_3 = {"fileSave", 0, 0 };
    static const QUParameter param_slot_4[] = {
	{ "filename", &static_QUType_QString, 0, QUParameter::In }
    };
    static const QUMethod slot_4 = {"fileSaveAs", 1, param_slot_4 };
    static const QUMethod slot_5 = {"filePrint", 0, 0 };
    static const QUMethod slot_6 = {"fileExit", 0, 0 };
    static const QUMethod slot_7 = {"helpIndex", 0, 0 };
    static const QUMethod slot_8 = {"helpContents", 0, 0 };
    static const QUMethod slot_9 = {"helpAbout", 0, 0 };
    static const QUMethod slot_10 = {"constructDual", 0, 0 };
    static const QUMethod slot_11 = {"updateMesh", 0, 0 };
    static const QUMethod slot_12 = {"DebugButton_pressed", 0, 0 };
    static const QUMethod slot_13 = {"ActorTreeView_selectionChanged", 0, 0 };
    static const QUMethod slot_14 = {"TagTreeView_selectionChanged", 0, 0 };
    static const QUMethod slot_15 = {"displayVisible", 0, 0 };
    static const QUMethod slot_16 = {"displayDraw", 0, 0 };
    static const QUMethod slot_17 = {"displayWireframeShaded", 0, 0 };
    static const QUMethod slot_18 = {"displayInvertSelection", 0, 0 };
    static const QUParameter param_slot_19[] = {
	{ "item", &static_QUType_ptr, "QListViewItem", QUParameter::In },
	{ 0, &static_QUType_varptr, "\x0e", QUParameter::In },
	{ 0, &static_QUType_int, 0, QUParameter::In }
    };
    static const QUMethod slot_19 = {"ActorListView1_rightButtonPressed", 3, param_slot_19 };
    static const QUMethod slot_20 = {"displayInvisible", 0, 0 };
    static const QUParameter param_slot_21[] = {
	{ 0, &static_QUType_ptr, "QListViewItem", QUParameter::In },
	{ 0, &static_QUType_varptr, "\x0e", QUParameter::In },
	{ 0, &static_QUType_int, 0, QUParameter::In }
    };
    static const QUMethod slot_21 = {"TagListView1_rightButtonPressed", 3, param_slot_21 };
    static const QUMethod slot_22 = {"displayDisplayAll", 0, 0 };
    static const QUMethod slot_23 = {"CropToolButton_clicked", 0, 0 };
    static const QUMethod slot_24 = {"displayDrawSheetAction_activated", 0, 0 };
    static const QUMethod slot_25 = {"APbutton_clicked", 0, 0 };
    static const QUMethod slot_26 = {"negAPbutton_clicked", 0, 0 };
    static const QUMethod slot_27 = {"FOCbutton_clicked", 0, 0 };
    static const QUMethod slot_28 = {"FSbutton_clicked", 0, 0 };
    static const QUMethod slot_29 = {"negFCbutton_clicked", 0, 0 };
    static const QUMethod slot_30 = {"fileSaveAs", 0, 0 };
    static const QUMethod slot_31 = {"resetDisplay", 0, 0 };
    static const QUMethod slot_32 = {"redrawDisplay", 0, 0 };
    static const QUMethod slot_33 = {"languageChange", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "fileNew()", &slot_0, QMetaData::Public },
	{ "fileOpen()", &slot_1, QMetaData::Public },
	{ "fileOpen(const QString&)", &slot_2, QMetaData::Public },
	{ "fileSave()", &slot_3, QMetaData::Public },
	{ "fileSaveAs(const QString&)", &slot_4, QMetaData::Public },
	{ "filePrint()", &slot_5, QMetaData::Public },
	{ "fileExit()", &slot_6, QMetaData::Public },
	{ "helpIndex()", &slot_7, QMetaData::Public },
	{ "helpContents()", &slot_8, QMetaData::Public },
	{ "helpAbout()", &slot_9, QMetaData::Public },
	{ "constructDual()", &slot_10, QMetaData::Public },
	{ "updateMesh()", &slot_11, QMetaData::Public },
	{ "DebugButton_pressed()", &slot_12, QMetaData::Public },
	{ "ActorTreeView_selectionChanged()", &slot_13, QMetaData::Public },
	{ "TagTreeView_selectionChanged()", &slot_14, QMetaData::Public },
	{ "displayVisible()", &slot_15, QMetaData::Public },
	{ "displayDraw()", &slot_16, QMetaData::Public },
	{ "displayWireframeShaded()", &slot_17, QMetaData::Public },
	{ "displayInvertSelection()", &slot_18, QMetaData::Public },
	{ "ActorListView1_rightButtonPressed(QListViewItem*,const QPoint&,int)", &slot_19, QMetaData::Public },
	{ "displayInvisible()", &slot_20, QMetaData::Public },
	{ "TagListView1_rightButtonPressed(QListViewItem*,const QPoint&,int)", &slot_21, QMetaData::Public },
	{ "displayDisplayAll()", &slot_22, QMetaData::Public },
	{ "CropToolButton_clicked()", &slot_23, QMetaData::Public },
	{ "displayDrawSheetAction_activated()", &slot_24, QMetaData::Public },
	{ "APbutton_clicked()", &slot_25, QMetaData::Public },
	{ "negAPbutton_clicked()", &slot_26, QMetaData::Public },
	{ "FOCbutton_clicked()", &slot_27, QMetaData::Public },
	{ "FSbutton_clicked()", &slot_28, QMetaData::Public },
	{ "negFCbutton_clicked()", &slot_29, QMetaData::Public },
	{ "fileSaveAs()", &slot_30, QMetaData::Public },
	{ "resetDisplay()", &slot_31, QMetaData::Public },
	{ "redrawDisplay()", &slot_32, QMetaData::Public },
	{ "languageChange()", &slot_33, QMetaData::Protected }
    };
    static const QUMethod signal_0 = {"toggled", 0, 0 };
    static const QMetaData signal_tbl[] = {
	{ "toggled()", &signal_0, QMetaData::Public }
    };
    metaObj = QMetaObject::new_metaobject(
	"uiQVDual", parentObject,
	slot_tbl, 34,
	signal_tbl, 1,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_uiQVDual.setMetaObject( metaObj );
    return metaObj;
}

void* uiQVDual::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "uiQVDual" ) )
	return this;
    return QMainWindow::qt_cast( clname );
}

// SIGNAL toggled
void uiQVDual::toggled()
{
    activate_signal( staticMetaObject()->signalOffset() + 0 );
}

bool uiQVDual::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: fileNew(); break;
    case 1: fileOpen(); break;
    case 2: fileOpen((const QString&)static_QUType_QString.get(_o+1)); break;
    case 3: fileSave(); break;
    case 4: fileSaveAs((const QString&)static_QUType_QString.get(_o+1)); break;
    case 5: filePrint(); break;
    case 6: fileExit(); break;
    case 7: helpIndex(); break;
    case 8: helpContents(); break;
    case 9: helpAbout(); break;
    case 10: constructDual(); break;
    case 11: updateMesh(); break;
    case 12: DebugButton_pressed(); break;
    case 13: ActorTreeView_selectionChanged(); break;
    case 14: TagTreeView_selectionChanged(); break;
    case 15: displayVisible(); break;
    case 16: displayDraw(); break;
    case 17: displayWireframeShaded(); break;
    case 18: displayInvertSelection(); break;
    case 19: ActorListView1_rightButtonPressed((QListViewItem*)static_QUType_ptr.get(_o+1),(const QPoint&)*((const QPoint*)static_QUType_ptr.get(_o+2)),(int)static_QUType_int.get(_o+3)); break;
    case 20: displayInvisible(); break;
    case 21: TagListView1_rightButtonPressed((QListViewItem*)static_QUType_ptr.get(_o+1),(const QPoint&)*((const QPoint*)static_QUType_ptr.get(_o+2)),(int)static_QUType_int.get(_o+3)); break;
    case 22: displayDisplayAll(); break;
    case 23: CropToolButton_clicked(); break;
    case 24: displayDrawSheetAction_activated(); break;
    case 25: APbutton_clicked(); break;
    case 26: negAPbutton_clicked(); break;
    case 27: FOCbutton_clicked(); break;
    case 28: FSbutton_clicked(); break;
    case 29: negFCbutton_clicked(); break;
    case 30: fileSaveAs(); break;
    case 31: resetDisplay(); break;
    case 32: redrawDisplay(); break;
    case 33: languageChange(); break;
    default:
	return QMainWindow::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool uiQVDual::qt_emit( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->signalOffset() ) {
    case 0: toggled(); break;
    default:
	return QMainWindow::qt_emit(_id,_o);
    }
    return TRUE;
}
#ifndef QT_NO_PROPERTIES

bool uiQVDual::qt_property( int id, int f, QVariant* v)
{
    return QMainWindow::qt_property( id, f, v);
}

bool uiQVDual::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
