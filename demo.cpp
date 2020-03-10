/*
  ======================================================================
   demo.c --- protoype to show off the simple solver
  ----------------------------------------------------------------------
   Author : Jos Stam (jstam@aw.sgi.com)
   Creation Date : Jan 9 2003

   Description:

	This code is a simple prototype that demonstrates how to use the
	code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
	for Games". This code uses OpenGL and GLUT for graphics and interface

  =======================================================================
*/

#include <stdlib.h>
#include <stdio.h>
#include "project.h"

/* global variables */

static float dt, diff, visc;
static float force, source;
static uint32_t dvel;

static float * u, * v, * u_prev, * v_prev;
static float * dens, * dens_prev;

static uint32_t win_id;
static int win_x, win_y;
static uint32_t mouse_down[3];
static int omx, omy, mx, my;

sf::RenderWindow *window;
sf::Clock *g_clock;
sf::Font *font;

float timeSpeed;
const uint32_t frameRate = 60;
const float visLen = 30.0f;
const uint32_t resolution = 512;
uint32_t iterations;
uint32_t N;

bool profiling;
int optim_mode;

/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/

static void free_data ( void )
{
	if ( u ) free ( u );
	if ( v ) free ( v );
	if ( u_prev ) free ( u_prev );
	if ( v_prev ) free ( v_prev );
	if ( dens ) free ( dens );
	if ( dens_prev ) free ( dens_prev );
}

static void clear_data ( void )
{
	uint32_t i, size=(N+2)*(N+2);

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
	}
}

static uint32_t allocate_data ( void )
{
	uint32_t size = (N+2)*(N+2);

	u			= (float *) malloc ( size*sizeof(float) );
	v			= (float *) malloc ( size*sizeof(float) );
	u_prev		= (float *) malloc ( size*sizeof(float) );
	v_prev		= (float *) malloc ( size*sizeof(float) );
	dens		= (float *) malloc ( size*sizeof(float) );	
	dens_prev	= (float *) malloc ( size*sizeof(float) );

	if ( !u || !v || !u_prev || !v_prev || !dens || !dens_prev ) {
		fprintf ( stderr, "cannot allocate data\n" );
		return ( 0 );
	}

	return ( 1 );
}


/*
  ----------------------------------------------------------------------
   SFML specific routines
  ----------------------------------------------------------------------
*/

static void init_sfml(uint32_t resolution)
{	
	font = new sf::Font();
	font->loadFromFile("UbuntuMono-B.ttf");
	sf::ContextSettings settings;
	settings.antialiasingLevel = 4;
	window = new sf::RenderWindow(sf::VideoMode(resolution, resolution), "Wind Field Prototype", sf::Style::Default, settings);
	window->setFramerateLimit(frameRate);
	g_clock = new sf::Clock();
}

static void mouse_func(uint32_t id, bool pressed)
{
	mouse_down[id] = pressed;
	
	sf::Vector2i mPos = sf::Mouse::getPosition(*window);

	omx = mx = mPos.x;
	omx = my = mPos.y;

}

static void motion_func()
{
	sf::Vector2i mPos = sf::Mouse::getPosition(*window);
	mx = mPos.x;
	my = mPos.y;
}

static void draw_density ( void )
{
	// Density rendering not yet implemented
}


static void get_from_UI(float *d, float *u, float *v)
{
	uint32_t i, j, size = (N + 2) * (N + 2);

	for (i = 0; i < size; i++)
	{
		u[i] = v[i] = d[i] = 0.0f;
	}

	if (!mouse_down[0] && !mouse_down[2])
		return;

	i = (int)((mx / (float)win_x) * N + 1);
	j = (int)(((win_y - my) / (float)win_y) * N + 1);

	if (i < 1 || i > N || j < 1 || j > N)
		return;

	if (mouse_down[0]){
		if(optim_mode){
			u[ZIX(i, j)] = force * dt * (mx - omx);
			v[ZIX(i, j)] = force * dt * (omy - my);
		}
		else{
			u[IX(i, j)] = force * dt * (mx - omx);
			v[IX(i, j)] = force * dt * (omy - my);
		}
	}

	if (mouse_down[2]){
		if (optim_mode){
			d[ZIX(i, j)] = source;
		}
		else{
			d[IX(i, j)] = source;
		}
	}

	omx = mx;
	omy = my;

	return;
}

static void ProfileLoop()
{
	while(iterations > 0){
		iterations--;
		switch(optim_mode){
			case 0:
			{
				// Constant force is added each iteration to demonstrate simulation without needing input.
				base::add_force(N / 4, N / 4, 100.0f, 0); 
				base::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				base::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				base::dens_step(N, dens, dens_prev, u, v, diff, dt);
				break;
			}
			case 1:
			{
				// Constant force is added each iteration to demonstrate simulation without needing input.
				opt::add_force(N / 4, N / 4, 100.0f, 0); 
				opt::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				opt::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				opt::dens_step(N, dens, dens_prev, u, v, diff, dt);
				break;
			}
			case 2:
			{
				// Constant force is added each iteration to demonstrate simulation without needing input.
				opt::add_force(N / 4, N / 4, 100.0f, 0); 
				opt::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				parallel::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				parallel::dens_step(N, dens, dens_prev, u, v, diff, dt);
				break;
			}
		}
	}

}

static void GameLoop(){
	DPRINT("GameLoop begin\n");
	while (window->isOpen())
	{
		sf::Event event;
		while (window->pollEvent(event))
		{
			if (event.type == sf::Event::Closed){
				DPRINT("Closed\n");
				window->close();
			}
			if (event.type == sf::Event::MouseButtonPressed){
				if(event.mouseButton.button == sf::Mouse::Left){
					mouse_func(0, true);
				}
				if (event.mouseButton.button == sf::Mouse::Middle){
					mouse_func(1, true);
				}
				if (event.mouseButton.button == sf::Mouse::Right){
					mouse_func(2, true);
				}
			}
			if (event.type == sf::Event::MouseButtonReleased){
				if (event.mouseButton.button == sf::Mouse::Left){
					mouse_func(0, false);
				}
				if (event.mouseButton.button == sf::Mouse::Middle){
					mouse_func(1, false);
				}
				if (event.mouseButton.button == sf::Mouse::Right){
					mouse_func(2, false);
				}
			}
		}

		if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
		{
			motion_func();
		}

		get_from_UI(dens_prev, u_prev, v_prev);
	
		DPRINT("Iteration "<< iterations <<"\n");
		switch(optim_mode){
			case 0:
			{
				base::add_force(N / 4, N / 4, 100.0f, 0); // Constant force is added each iteration to demonstrate simulation without needing input.
				base::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				base::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				base::dens_step(N, dens, dens_prev, u, v, diff, dt);
				base::render_velocity();
				break;
			}
			case 1:
			{
				opt::add_force(N / 4, N / 4, 100.0f, 0); // Constant force is added each iteration to demonstrate simulation without needing input.
				opt::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				opt::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				opt::dens_step(N, dens, dens_prev, u, v, diff, dt);
				opt::render_velocity();
				break;
			}
			case 2:
			{
				opt::add_force(N / 4, N / 4, 100.0f, 0); // Constant force is added each iteration to demonstrate simulation without needing input.
				opt::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				parallel::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				parallel::dens_step(N, dens, dens_prev, u, v, diff, dt);
				opt::render_velocity();
				break;
			}
		}
	}
}

namespace opt
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		if (i < 1 || i > N || j < 1 || j > N)
			return;
		u[ZIX(i, j)] += dt * xForce;
		v[ZIX(i, j)] += dt * yForce;
	}
	void render_velocity(){
		DPRINT("render_velocity begin\n");
		sf::Vertex line[2*(N+2)*(N+2)];
		DPRINT("1\n");
		#ifdef DISPLAY_COORDS
		sf::Text texts[(N+2)*(N+2)];
		#endif
		DPRINT("2\n");
		bool hit[2*(N+2)*(N+2)];
		DPRINT("3\n");
		uint32_t hits = 0;
		uint32_t zX, zY, i, j, k;

		DPRINT("RENDERLOOP: zonesInRow = (" << zonesInRow <<")\n");
		DPRINT("RENDERLOOP: zoneLen = (" << zoneLen <<")\n");
		DPRINT("RENDERLOOP: ARRAY SIZE = (" << 2*N*N <<")\n");


		for( zY=0; zY<zonesInRow; zY++){
			for( zX=0; zX<zonesInRow; zX++){
				for ( j=zY*zoneLen; j<(zY+1)*zoneLen; j++ ) {
					for ( i=zX*zoneLen; i<(zX+1)*zoneLen; i++ ) { 
						if(i <= 0 || j <= 0 || i >= N || j >= N){
							DPRINT("RENDERLOOP: (" << i << "," << j <<") OMITTED\n");	
							#ifdef DISPLAY_COORDS				
							uint32_t zix = ZIX(i, j);
							float xStep = (float)resolution/N;
							float yStep = (float)resolution/N;	
							texts[zix] = sf::Text(sf::String("("+std::to_string(i)+","+std::to_string(j)+")"), *font);
							texts[zix].setCharacterSize(11);
							texts[zix].setFillColor(sf::Color::Red);
							sf::FloatRect textRect = texts[zix].getLocalBounds();
							float xPos = xStep * ((float)i);
							float yPos = resolution-(yStep * ((float)j));
							xPos -= textRect.width/2;
							yPos -= textRect.height/2;
							texts[zix].setPosition(xPos, yPos);
							#endif
							continue;
						}

						uint32_t zix = ZIX(i, j);
						uint32_t index = 2 * zix;
						DPRINT("RENDERLOOP: (" << i << "," << j <<")="<< index <<"\n");
						
						hit[index] = true;
						hit[index+1] = true;
						hits++;

						float xStep = (float)resolution/N;
						float yStep = (float)resolution/N;

						float xPos = xStep * ((float)i);
						float yPos = resolution-(yStep * ((float)j));
						float xPerc = (float)i / (float)N;
						float yPerc = (float)j / (float)N;
						float xOffset = u[zix]*visLen;
						float yOffset = v[zix]*visLen;

						float colInterp = (abs(xOffset)+abs(yOffset))/(visLen*0.5f);
						colInterp = colInterp <= 1 ? colInterp : 1.0f;

						line[index].position = sf::Vector2f(xPos, yPos);
						line[index].color.r = 255*colInterp;
						line[index].color.b = 255*(1-colInterp);
						line[index].color.g = 128;
						line[index+1].position = sf::Vector2f(xPos + xOffset,  yPos + yOffset);
						line[index+1].color = line[index].color;
						#ifdef DISPLAY_COORDS
						texts[zix] = sf::Text(sf::String("("+std::to_string(i)+","+std::to_string(j)+")"), *font);
						texts[zix].setCharacterSize(11);
						texts[zix].setFillColor(line[index].color);
						sf::FloatRect textRect = texts[zix].getLocalBounds();
						texts[zix].setPosition(xPos - textRect.width/2, yPos- textRect.height/2);
						#endif
					}
				}	
			}
		}

		for(uint32_t i = 0; i < 2 * N * N; i++)
		{
			if(hit[i] == false){
				DPRINT("Cell at "<<i<<" was missed by new indexing!\n");
			}
		}

		window->clear();
		DPRINT("Cleared window \n");
		#ifdef DISPLAY_COORDS
		if(DISPLAY_COORDS){
			for(uint32_t i = 0; i < (N+2)*(N+2); i++){
				window->draw(texts[i]);
			}
		}
		#endif
		window->draw(line, 2 * (N+2) * (N+2), sf::PrimitiveType::Lines);
		window->display();
	}
} // namespace opt

namespace base
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		if (i < 1 || i > N || j < 1 || j > N)
			return;
		u[IX(i, j)] += dt * xForce;
		v[IX(i, j)] += dt * yForce;
	}
	void render_velocity(){
		sf::Vertex line[2*N*N];
		for(uint32_t y = 0; y < N; y++){
			for(uint32_t x = 0; x < N; x++){
				uint32_t index = 2*(x+(y*N));
				//DPRINT("("<< x*2 << ", " << y << ")=ArrayPos " << index << "\n"); 
				//DPRINT("("<< x*2+1 << ", " << y << ")=ArrayPos " << index+1 << "\n");
				
				float xStep = (float)resolution/N;
				float yStep = (float)resolution/N;

				float xPos = xStep * ((float)x + 0.5f);
				float yPos = yStep * ((float)y + 0.5f);
				float xPerc = (float)x / (float)N;
				float yPerc = (float)y / (float)N;				

				float xOffset = u[IX(x, N-y)]*visLen;
				float yOffset = v[IX(x, N-y)]*visLen;
				float colInterp = (abs(xOffset)+abs(yOffset))/(visLen*0.5f);
				colInterp = colInterp <= 1 ? colInterp : 1.0f;

				line[index].color.r = 255*colInterp;
				line[index].color.b = 255*(1-colInterp);
				line[index].color.g = 128;
				
				line[index].position = sf::Vector2f(xPos, yPos);

				line[index+1].position = sf::Vector2f(xPos+xOffset, yPos+yOffset);
				line[index+1].color = line[index].color;
			}
		}
		window->clear();
		window->draw(line, 2 * N * N, sf::PrimitiveType::Lines);
		window->display();
	}
} // namespace base

/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{

	if ( argc != 1 && argc != 5 ) {
		fprintf ( stderr, "Usage: demo.exe <profiling mode[0,1]> <use optimizations[0,1]> <resolution[int]> <iterations[int]>\n");
		return 1;
	}

	if ( argc == 1 ) {
		profiling = false;
		optim_mode = false;
		N = 126;
		iterations = 600;
		fprintf ( stderr, "Using defaults: visualizer mode, unoptim_mode, 128x128\n");
	} else {
		profiling = atoi(argv[1]) > 0 ? true : false;
		optim_mode = atoi(argv[2]);
		N = atoi(argv[3])-2;
		iterations = atoi(argv[4]);
	}
	timeSpeed = 1.0f;
	diff = 0.0f;
	visc = 0.0f;
	force = 1000.0f;
	source = 100.0f;

	//zoneSize = zoneLen*zoneLen;
	//zoneDivisor = 1/(float)zoneLen;
	//zonesInRow = N/zoneLen;

	DPRINT("zoneSize:" << zoneSize <<"\n");
	//DPRINT("zoneDivisor:" << zoneDivisor <<"\n");
	DPRINT("zonesInRow:" << zonesInRow <<"\n");

	printf ( "\n\nHow to use this demo:\n\n" );
	//printf ( "\t Add densities with the right mouse button\n" );
	printf ( "Add velocities with the left mouse button and dragging the mouse\n" );
	//printf ( "\t Toggle density/velocity display with the 'v' key\n" );
	//printf ( "\t Clear the simulation by pressing the 'c' key\n" );
	//printf ( "\t Quit by pressing the 'q' key\n" );

	dvel = 0;

	if ( !allocate_data () ) return 1;
	clear_data ();

	win_x = resolution;
	win_y = resolution;

	DPRINT("Max index of ZIX: " << ZIX(N, N) << "\nMax size of array:" << (N+2)*(N+2) << ".\n");

	if (profiling)
	{
		dt = timeSpeed / 60.0f;
	}
	else{
		dt = timeSpeed / (float)frameRate;
	}
	if(profiling){
		ProfileLoop();
	}
	else{
		init_sfml(resolution);
		GameLoop();
	}
	return 0;
}