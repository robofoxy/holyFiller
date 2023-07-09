#include <iostream>
#include <queue>
#include "Painter.h"
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

int main()
{
	HWND window = SoWin::init("");

	SoWinExaminerViewer* viewer = new SoWinExaminerViewer(window);

	//make a dead simple scene graph by using the Coin library, only containing a single cone under the scenegraph root
	SoSeparator* root = new SoSeparator;
	root->ref();

	Painter* painter = new Painter();
	Mesh* mesh = new Mesh();

	std::string name = "bunny-1.obj";
	mesh->loadObj(name.c_str());

	vector<Eigen::Vector3i> filled;
	filled = holyFiller(mesh, ANGLE, UNIFORM);
	
	SoSeparator* thickEdgeSep = new SoSeparator();

	SoSeparator* sep = painter->getShapeSep(mesh, filled, thickEdgeSep);

	root->addChild(sep);
	root->addChild(thickEdgeSep);

	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	char k;
	scanf_s(" %c", &k);
	
	return 0;
}