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
/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you wish to add, delete or rename functions or slots use
** Qt Designer which will update this file, preserving your code. Create an
** init() function in place of a constructor, and a destroy() function in
** place of a destructor.
*****************************************************************************/


void uiQVDual::fileNew()
{
  
}


void uiQVDual::fileOpen()
{
  assert(NULL != vtkMOABUtils::mbImpl);
  QString filename = QFileDialog::getOpenFileName(
      QString::null,
      "Exodus files (*.g *.gen *.exoII);;Cub files (*.cub);;SLAC Netcdf files (*.ncdf);;All Files (*)", this
      );
  if(filename == QString::null)
    return;

  fileOpen(filename);
}

void uiQVDual::fileOpen( const QString &filename )
{
  // create a MOAB vtk reader and pass responsibility to that
  vtkMOABReader *reader = vtkMOABReader::New();
  
  // if we already have a grid, make sure it's that reader's output too
  if (NULL != vtkMOABUtils::myUG) reader->SetOutput(vtkMOABUtils::myUG);
  
  reader->SetFileName(filename.ascii());
  
  // need to update here, in case we're doing something else which requires moab data
  reader->Update();
  
  // keep track of the reader's output, to know where the cell and point data is
  vtkMOABUtils::myUG = reader->GetOutput();

    // update anything else in the UI
  this->updateMesh();
  
    // Reset camera
  vtkMOABUtils::myRen->ResetCamera();
}

void uiQVDual::fileSave()
{
  
}


void uiQVDual::fileSaveAs()
{
  
}


void uiQVDual::filePrint()
{
  
}


void uiQVDual::fileExit()
{
  qApp->exit();
}


void uiQVDual::helpIndex()
{
  
}


void uiQVDual::helpContents()
{
  
}


void uiQVDual::helpAbout()
{
  
}


void uiQVDual::init()
{
  // initialize vtk stuff
  vtkRenderer *ren = vtkRenderer::New();
  vtkWidget->GetRenderWindow()->AddRenderer(ren);
  
  // initialize MOAB if it's not already
  if (NULL == vtkMOABUtils::mbImpl) 
    vtkMOABUtils::init(NULL, ren);
  
  currentWin = 0;

  cropToolPopup = NULL;

  drawDual = NULL;
}


void uiQVDual::destroy()
{
  vtkMOABUtils::destroy();
}



void uiQVDual::constructDual()
{
  // tell MOAB to construct the dual first
  DualTool dt(vtkMOABUtils::mbImpl);
  MBErrorCode result = dt.construct_hex_dual();
  if (MB_SUCCESS == result) {
    // success - now populate vtk data; first the points
    result = vtkMOABUtils::make_vertex_points(vtkMOABUtils::myUG);
    if (MB_SUCCESS != result) return;
    
    // now the polylines
    result = vtkMOABUtils::make_cells(MBEDGE, vtkMOABUtils::myUG);
    if (MB_SUCCESS != result) return;
    
    // now the polygons
    result = vtkMOABUtils::make_cells(MBPOLYGON, vtkMOABUtils::myUG);
    if (MB_SUCCESS != result) return;
    
    // finally, the sets
    MBRange dual_sets;
    result = dt.get_dual_hyperplanes(vtkMOABUtils::mbImpl, 1, dual_sets);
    if (MB_SUCCESS != result) return;
    result = vtkMOABUtils::update_set_actors(dual_sets, vtkMOABUtils::myUG, true, true, true);
    if (MB_SUCCESS != result) return;

    dual_sets.clear();
    result = dt.get_dual_hyperplanes(vtkMOABUtils::mbImpl, 2, dual_sets);
    if (MB_SUCCESS != result) return;
    result = vtkMOABUtils::update_set_actors(dual_sets, vtkMOABUtils::myUG, true, false, true);

    vtkMOABUtils::construct_lookup_table(dual_sets.size());
  }
  
  updateMesh();
}

void uiQVDual::updateMesh()
{
  // mesh was updated; update the various UI elements
  this->updateTagList();
  this->updateActorList();
  
  // Render
  vtkWidget->GetRenderWindow()->Render();
}



void uiQVDual::DebugButton_pressed()
{
  vtkMOABUtils::print_debug();
}


void uiQVDual::updateTagList()
{
  // clean the tag list view
  TagListView1->clear();    
  
  // get a list of tags
  std::vector<MBTag> tag_handles;
  MBErrorCode result = vtkMOABUtils::mbImpl->tag_get_tags(tag_handles);
  if (MB_SUCCESS != result) return;

  // display each tag as a parent item
  QListViewItem *tags_item;
  
  for (std::vector<MBTag>::iterator tag_it = tag_handles.begin(); tag_it != tag_handles.end();
       tag_it++) {
    
    std::string tag_name;
    result = vtkMOABUtils::mbImpl->tag_get_name(*tag_it, tag_name);
    if (MB_SUCCESS != result) continue;

    // don't display tags with "__" prefix
    if (0 == strncmp(tag_name.c_str(), "__", 2))
      continue;
    
    // get all the sets which contain this tag
    MBRange tag_sets;
    result = vtkMOABUtils::mbImpl->get_entities_by_type_and_tag(
      0, MBENTITYSET, &(*tag_it), NULL, 1, tag_sets, MBInterface::UNION);
    if (MB_SUCCESS != result) continue;

    // create a parent tag item
    tags_item = new QListViewItem(TagListView1, tag_name.c_str());
    tags_item->setOpen(false);
    itemSetMap[tags_item] = 0;

    int i;
    MBRange::iterator set_it;
    QString set_id;
    QListViewItem *set_item;
    for (i = 0, set_it = tag_sets.begin(); set_it != tag_sets.end(); set_it++, i++) {
    // make an item for this set
      char set_name[CATEGORY_TAG_NAME_LENGTH];
      result = vtkMOABUtils::get_set_category_name(*set_it, set_name);
      if (MB_SUCCESS != result) sprintf(set_name, "(none)\0");
      set_item = new QListViewItem(tags_item, set_name);
      itemSetMap[set_item] = *set_it;
    }
  }
}


void uiQVDual::updateActorList()
{
    // clean the ActorView
  ActorListView1->clear();    

    // update "contains" view
  MBRange top_sets;
  QListViewItem *last_item;
  MBErrorCode result = vtkMOABUtils::get_top_contains_sets(top_sets);
  if (MB_SUCCESS == result && !top_sets.empty()) {
  
      // keep track of last actor line shown
    last_item = new QListViewItem(ActorListView1, "Contains Sets");
    itemSetMap[last_item] = 0;
    last_item->setOpen(false);

    for (MBRange::iterator rit = top_sets.begin(); rit != top_sets.end(); rit++)
      updateActorContainsList(last_item, *rit);
  }
  
    // update "parent/child" view
  top_sets.clear();
  result = vtkMOABUtils::get_top_parent_sets(top_sets);
  if (MB_SUCCESS == result && !top_sets.empty()) {
  
      // keep track of last actor line shown
    last_item = new QListViewItem(ActorListView1, "Parent/Child Sets");
    itemSetMap[last_item] = 0;
    last_item->setOpen(false);

    for (MBRange::iterator rit = top_sets.begin(); rit != top_sets.end(); rit++)
      updateActorParentList(last_item, *rit);
  }
}



void uiQVDual::ActorTreeView_selectionChanged()
{
  currentWin = 1;

    // go through selected items, acting accordingly
  static std::set<QListViewItem*> high_items, unhigh_items, tmp_sel, tmp_unsel;
  high_items.clear(); unhigh_items.clear(); tmp_sel.clear(); tmp_unsel.clear();
  getSelected(ActorListView1, tmp_sel, tmp_unsel);
  // newly selected is difference between selected and itemSelList
  std::set_difference(tmp_sel.begin(), tmp_sel.end(),
                      itemSelList.begin(), itemSelList.end(),
                      std::inserter(high_items, high_items.begin()));
  
  // newly unselected is intersection between unselected and itemSelList
  std::set_intersection(tmp_unsel.begin(), tmp_unsel.end(),
                        itemSelList.begin(), itemSelList.end(),
                        std::inserter(unhigh_items, unhigh_items.begin()));

  if (!high_items.empty() || !unhigh_items.empty()) {
      // reset the properties of the items
    changeSetProperty(high_items, unhigh_items);
    
      // Re-render
    vtkWidget->GetRenderWindow()->Render();
  }
}

void uiQVDual::TagTreeView_selectionChanged()
{
  currentWin = 2;
  
    // go through selected items, acting accordingly
  static std::set<QListViewItem*> high_items, unhigh_items, tmp_sel, tmp_unsel;
  high_items.clear(); unhigh_items.clear(); tmp_sel.clear(); tmp_unsel.clear();
  getSelected(TagListView1, tmp_sel, tmp_unsel);
  // newly selected is difference between selected and itemSelList
  std::set_difference(tmp_sel.begin(), tmp_sel.end(),
                      itemSelList.begin(), itemSelList.end(),
                      std::inserter(high_items, high_items.begin()));
  
  // newly unselected is intersection between unselected and itemSelList
  std::set_intersection(tmp_unsel.begin(), tmp_unsel.end(),
                        itemSelList.begin(), itemSelList.end(),
                        std::inserter(unhigh_items, unhigh_items.begin()));
  
  if (!high_items.empty() || !unhigh_items.empty()) {
      // reset the properties of the items
    changeSetProperty(high_items, unhigh_items);
    
      // Re-render
    vtkWidget->GetRenderWindow()->Render();
  }
}

void uiQVDual::updateActorContainsList( QListViewItem *item, MBEntityHandle set_handle )
{
  vtkActor *set_actor;
  QListViewItem *set_item;

    // get the actor for this set
  MBErrorCode result = vtkMOABUtils::mbImpl->tag_get_data(vtkMOABUtils::vtkSetActorTag,
                                                          &set_handle, 1, &set_actor);
  if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) return;
  
  int num_sets;
  if (NULL == set_actor) {
      // no entities - check for contained sets, and if none, return
    result = vtkMOABUtils::mbImpl->get_number_entities_by_type(set_handle, MBENTITYSET,
                                                               num_sets);
    if (MB_SUCCESS != result || 0 == num_sets) return;
  }
  
    // has an actor, or contains sets; allocate a list item
  char set_name[CATEGORY_TAG_NAME_LENGTH];
  vtkMOABUtils::get_set_category_name(set_handle, set_name);
  set_item = new QListViewItem(item, set_name);
  itemSetMap[set_item] = set_handle;
  set_item->setOpen(false);

  if (NULL == set_actor && 0 == num_sets) return;
  
    // get the list of sets in this set
  MBRange contained_sets;
  result = vtkMOABUtils::mbImpl->get_entities_by_type(set_handle, MBENTITYSET, 
                                                      contained_sets);
  if (MB_SUCCESS != result) return;

  for (MBRange::iterator rit = contained_sets.begin(); rit != contained_sets.end(); rit++)
    updateActorContainsList(set_item, *rit);
}


void uiQVDual::updateActorParentList( QListViewItem *item, MBEntityHandle set_handle )
{
  vtkActor *set_actor;
  QListViewItem *set_item;

    // get the actor for this set
  MBErrorCode result = vtkMOABUtils::mbImpl->tag_get_data(vtkMOABUtils::vtkSetActorTag,
                                                          &set_handle, 1, &set_actor);
  if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) return;
  
  int num_sets;
  if (NULL == set_actor) {
      // no entities - check for child sets, and if none, return
    result = vtkMOABUtils::mbImpl->num_child_meshsets(set_handle, &num_sets);
    if (MB_SUCCESS != result || 0 == num_sets) return;
  }
  
    // has an actor, or parent sets; allocate a list item
  char set_name[CATEGORY_TAG_NAME_LENGTH];
  vtkMOABUtils::get_set_category_name(set_handle, set_name);
  set_item = new QListViewItem(item, set_name);
  itemSetMap[set_item] = set_handle;
  set_item->setOpen(false);

    // get the list of child sets
  if (NULL == set_actor && 0 == num_sets) return;

  std::vector<MBEntityHandle> child_sets;
  result = vtkMOABUtils::mbImpl->get_child_meshsets(set_handle, child_sets);
  if (MB_SUCCESS != result) return;

  for (std::vector<MBEntityHandle>::iterator vit = child_sets.begin(); 
       vit != child_sets.end(); vit++)
    updateActorParentList(set_item, *vit);
}

void uiQVDual::changeSetProperty( std::set<QListViewItem *> &high_sets, std::set<QListViewItem *> &unhigh_sets )
{
  MBRange high_mbsets, unhigh_mbsets;
  MBEntityHandle item_set;
  vtkActor *this_actor;
  
  for (std::set<QListViewItem *>::iterator sit = high_sets.begin(); 
       sit != high_sets.end(); sit++) {

    item_set = itemSetMap[*sit];
    if (0 == item_set)
        // must have child items - evaluate them too
      evalItem(*sit, true, high_mbsets, unhigh_mbsets);

    else {
      this_actor = vtkMOABUtils::get_actor(item_set);
      if (NULL == this_actor) evalItem(*sit, true, high_mbsets, unhigh_mbsets);
      else
        high_mbsets.insert(item_set);
    }
    
      // either way, add item to selection list
    itemSelList.insert(*sit);
  }
  for (std::set<QListViewItem *>::iterator sit = unhigh_sets.begin(); 
       sit != unhigh_sets.end(); sit++) {

    item_set = itemSetMap[*sit];
    if (0 == item_set)
        // must have child items - evaluate them too
      evalItem(*sit, false, unhigh_mbsets, unhigh_mbsets);

    else {
      this_actor = vtkMOABUtils::get_actor(item_set);
      if (NULL == this_actor) evalItem(*sit, true, high_mbsets, unhigh_mbsets);
      else
        unhigh_mbsets.insert(item_set);
    }

      // either way, remove item from selection list
    itemSelList.erase(*sit);
  }

  vtkMOABUtils::change_set_properties(high_mbsets, unhigh_mbsets);
}

void uiQVDual::evalItem( QListViewItem * item, const bool high, 
                         MBRange & high_mbsets, MBRange & unhigh_mbsets )
{
    // item is a parent, evaluate all the child items
  QListViewItem *child_item = item->firstChild();
  while (NULL != child_item) {
    MBEntityHandle child_set = itemSetMap[child_item];
    if (0 == child_set) 
        // if child item also has no set, must be a parent too
      this->evalItem(child_item, high, high_mbsets, unhigh_mbsets);
    else if (child_item->isSelected() != high) {
        // selection state is changing; put item set in proper (high or unhigh) set
      if (high) high_mbsets.insert(child_set);
      else unhigh_mbsets.insert(child_set);
        // since selection state changed on parent item and not child, need to set
        // selection state of child directly
      child_item->setSelected(high);
    }
    
    child_item = child_item->nextSibling();
  }
}


void uiQVDual::displayVisible()
{
  MBRange selected, unselected;
  getSelected(NULL, selected, unselected);

    // get child sets of selected ones, they should be drawn too
  std::vector<MBEntityHandle> children;
  for (MBRange::iterator rit = selected.begin(); rit != selected.end(); rit++)
    vtkMOABUtils::mbImpl->get_child_meshsets(*rit, children, 0);

  std::copy(children.begin(), children.end(), mb_range_inserter(selected));
  
  unselected.clear();
  vtkMOABUtils::change_set_visibility(selected, unselected);
  vtkWidget->GetRenderWindow()->Render();
}


void uiQVDual::displayDraw()
{
  MBRange selected, unselected;
  getSelected(NULL, selected, unselected);
  
    // get child sets of selected ones, they should be drawn too
  std::vector<MBEntityHandle> children;
  for (MBRange::iterator rit = selected.begin(); rit != selected.end(); rit++)
    vtkMOABUtils::mbImpl->get_child_meshsets(*rit, children, 0);

  std::copy(children.begin(), children.end(), mb_range_inserter(selected));
  
  vtkMOABUtils::change_set_visibility(selected, unselected);

    // unhighlight the visible ones; don't change the others
  unselected.clear();
  vtkMOABUtils::change_set_properties(unselected, selected);

  vtkWidget->GetRenderWindow()->Render();
}


void uiQVDual::displayWireframeShaded()
{
  MBRange selected, unselected;
  getSelected(NULL, selected, unselected);
  
  vtkMOABUtils::toggle_wireframe_shaded(selected);
}

void uiQVDual::displayInvertSelection()
{

}

void uiQVDual::ActorListView1_rightButtonPressed( QListViewItem *item, const QPoint &, int )
{
  currentWin = 1;
  Display->exec(QCursor::pos());
}

void uiQVDual::displayInvisible()
{
  MBRange selected, unselected;
  getSelected(NULL, selected, unselected);

    // get child sets of selected ones, they should be set too
  std::vector<MBEntityHandle> children;
  for (MBRange::iterator rit = selected.begin(); rit != selected.end(); rit++)
    vtkMOABUtils::mbImpl->get_child_meshsets(*rit, children, 0);

  std::copy(children.begin(), children.end(), mb_range_inserter(selected));
  
  unselected.clear();
  vtkMOABUtils::change_set_visibility(unselected, selected);
  vtkWidget->GetRenderWindow()->Render();
}

void uiQVDual::getSelected( QListView *listv, 
                            std::set<QListViewItem *> &selected, 
                            std::set<QListViewItem *> &unselected )
{
  if (NULL == listv && currentWin == 0) {
      // check all items
    for (std::map<QListViewItem*, MBEntityHandle>::iterator mit = itemSetMap.begin();
         mit != itemSetMap.end(); mit++) {
      if ((*mit).first->isSelected())
        selected.insert((*mit).first);
      else
        unselected.insert((*mit).first);
    }
    
    return;
  }
  
  QListViewItemIterator it;
  if (NULL == listv) {
    if (1 == currentWin) 
      it = QListViewItemIterator(ActorListView1);
    else if (2 == currentWin)
      it = QListViewItemIterator(TagListView1);
  }
  else it = QListViewItemIterator(listv);
  
  while ( it.current() ) {
    if ( it.current()->isSelected())
     selected.insert(*it);
    else
      unselected.insert(*it);
    ++it;
  }
}

void uiQVDual::getSelected( QListView *listv, MBRange &selected, MBRange &unselected )
{
  static std::set<QListViewItem*> selected_s, unselected_s;
  selected_s.clear(); unselected_s.clear();
  getSelected(listv, selected_s, unselected_s);
  getItemSets(selected_s, selected);
  getItemSets(unselected_s, unselected);
}

void uiQVDual::TagListView1_rightButtonPressed( QListViewItem *, const QPoint &, int )
{
  currentWin = 2;
  Display->exec(QCursor::pos());
}


void uiQVDual::displayDisplayAll()
{
  currentWin = 0;
  MBRange selected, unselected;
  getSelected(NULL, selected, unselected);
  selected.merge(unselected);
  selected.erase(1);
  unselected.clear();
  vtkMOABUtils::change_set_visibility(selected, unselected);
  vtkWidget->GetRenderWindow()->Render();

}

void uiQVDual::getItemSets( std::set<QListViewItem *> &items, MBRange &sets )
{
  std::set<QListViewItem*>::iterator sit;
  std::map<QListViewItem*, MBEntityHandle>::iterator mit;
  for (sit = items.begin(); sit != items.end(); sit++) {
    mit = itemSetMap.find(*sit);
    if (mit != itemSetMap.end() && (*mit).second != 0) 
      sets.insert((*mit).second);
  }
}


void uiQVDual::CropToolButton_clicked()
{
  if (NULL == cropToolPopup) {
    cropToolPopup = new CropToolPopup();
    cropToolPopup->vtk_widget(vtkWidget);
  }

  cropToolPopup->show();

}


void uiQVDual::displayDrawSheetAction_activated()
{
  MBRange selected, unselected;
  getSelected(NULL, selected, unselected);

    // get selected sets which are dual surfaces
  MBRange dual_surfs;
  DualTool dual_tool(vtkMOABUtils::mbImpl);
  MBEntityHandle dum_tag;
  MBErrorCode result;
  for (MBRange::iterator rit = selected.begin(); rit != selected.end(); rit++) {
    result = vtkMOABUtils::mbImpl->tag_get_data(dual_tool.dualSurface_tag(), &(*rit), 1,
                                                &dum_tag);
    if (MB_SUCCESS == result && 0 != dum_tag)
      dual_surfs.insert(*rit);
  }

    // now draw them
  if (NULL == drawDual) drawDual = new DrawDual();
  drawDual->draw_dual_surfs(dual_surfs);
}
