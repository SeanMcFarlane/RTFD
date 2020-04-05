/*
  ======================================================================
   demo.cpp - SFML-enabled version of the demo that provides visual feedback
  ----------------------------------------------------------------------
   Author : Sean Mcfarlane
   Creation Date : Jan 2019

   Description:

    This project was created based on the paper "Real-Time Fluid Dynamics for Games" by Jos Stam.
    Code is a modified version of the sample code provided alongside the paper, which was adapted to use C++ and SFML rather than C and OpenGL.
    This was created for the purpose of testing optimization methods as part of my CSC485C coursework.

  =======================================================================
*/

#include <stdlib.h>
#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "project.h"

/* global variables */

static uint32_t dvel;

static int win_x, win_y;
static uint32_t mouse_down[3];
static int omx, omy, mx, my;

sf::RenderWindow *window;
sf::Clock *g_clock;
sf::Font *font;

const float visLen = 30.0f;

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

// static void draw_density ( void )
// {
// 	// Density rendering not yet implemented
// }


static void Render(uint32_t optim_mode){
	switch(optim_mode){
		case 0: // Unoptimized
		{
			base::render_velocity();
			break;
		}
		case 1: // Single core optimized (ZIX)
		{
			opt::render_velocity();
			break;
		}
		case 2: // Parallelized (ZIX)
		{
			opt::render_velocity();
			break;
		}
		case 3: // SIMD (IX)
		{
			SIMD::render_velocity();
			break;
		}
		case 4: // SIMD & Parallelized (IX)
		{
			SIMD::render_velocity();
			break;
		}
	}
}

static void get_from_UI(float *d, float *u, float *v)
{
	uint32_t i, j, size = (N + bnd) * (N + bnd);

	for (i = 0; i < size; i++)
	{
		u[i] = v[i] = d[i] = 0.0f;
	}

	if (!mouse_down[0] && !mouse_down[2])
		return;

	i = (int)((mx / (float)win_x) * N)+pad;
	j = (int)(((win_y - my) / (float)win_y) * N)+pad;

	if (i < pad || i > N || j < pad || j > N)
		return;


	if (mouse_down[0]){
		if(optim_mode==1 || optim_mode==2){
			u[ZIX(i, j)] = force * dt * (mx - omx);
			v[ZIX(i, j)] = force * dt * (omy - my);
		}
		else{
			u[IX(i+pad, j+pad)] = force * dt * (mx - omx);
			v[IX(i+pad, j+pad)] = force * dt * (omy - my);
		}
	}

	if (mouse_down[2]){
		if(optim_mode==1 || optim_mode==2){
			d[ZIX(i, j)] = source;
		}
		else{
			d[IX(i+pad, j+pad)] = source;
		}
	}

	omx = mx;
	omy = my;

	return;
}

static void GameLoop(){
	DPRINT("GameLoop begin\n");
	while (window->isOpen())
	{
		cur_iter++;
        DPRINT("Iteration "<< cur_iter <<"\n");
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
		Simulate(optim_mode);
		Render(optim_mode);
	}
}

namespace opt
{
	void render_velocity(){
		DPRINT("render_velocity begin\n");
		sf::VertexArray line(sf::Lines, 2 * (N + 2) * (N + 2));

		DPRINT("1\n");
		#ifdef DISPLAY_COORDS
		sf::Text texts[(N+bnd)*(N+bnd)];
		#endif
		DPRINT("2\n");
		bool hit[2*(N+2)*(N+2)];
		DPRINT("3\n");
		uint32_t hits = 0;
		uint32_t zX, zY, i, j;

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
		window->draw(line);
		window->display();
	}
} // namespace opt

namespace base
{
	void render_velocity(){
		sf::VertexArray line(sf::Lines, 2 * N * N);
		for(uint32_t y = 0; y < N; y++){
			for(uint32_t x = 0; x < N; x++){
				uint32_t index = 2*(x+(y*N));
				//DPRINT("("<< x*2 << ", " << y << ")=ArrayPos " << index << "\n"); 
				//DPRINT("("<< x*2+1 << ", " << y << ")=ArrayPos " << index+1 << "\n");
				
				float xStep = (float)resolution/N;
				float yStep = (float)resolution/N;

				float xPos = xStep * ((float)x + 0.5f);
				float yPos = yStep * ((float)y + 0.5f);		

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
		window->draw(line);
		window->display();
	}
} // namespace base


namespace SIMD
{
	void render_velocity(){
		sf::VertexArray line(sf::Lines, 2 * N * N);
		for(uint32_t y = 0; y < N; y++){
			for(uint32_t x = 0; x < N; x++){
				uint32_t index = 2*(x+(y*N));
				//DPRINT("("<< x*2 << ", " << y << ")=ArrayPos " << index << "\n"); 
				//DPRINT("("<< x*2+1 << ", " << y << ")=ArrayPos " << index+1 << "\n");
				
				float xStep = (float)resolution/N;
				float yStep = (float)resolution/N;

				float xPos = xStep * ((float)x + 0.5f);
				float yPos = yStep * ((float)y + 0.5f);		

				float xOffset = u[IX(x+pad, N+pad-y)]*visLen;
				float yOffset = v[IX(x+pad, N+pad-y)]*visLen;
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
		window->draw(line);
		window->display();
	}
} // namespace SIMD


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{

	if ( argc != 1 && argc != 5 ) {
		fprintf ( stderr, "Usage: demo.exe <profiling mode[0-1]> <select implementation[0-4]> <resolution[int]> <iterations[int]>\n");
		return 1;
	}

	if ( argc == 1 ) {
		profiling = false;
		optim_mode = 4;
		N = 256;
		iterations = 1000;
		fprintf ( stderr, "Using defaults: visualizer mode, SIMD, 256x256\n");
	} else {
		profiling = atoi(argv[1]) > 0 ? true : false;
		optim_mode = atoi(argv[2]);
		N = atoi(argv[3]);
		iterations = atoi(argv[4]);
	}

	if(optim_mode==3||optim_mode==4){bnd = 8;}
	else{bnd = 2;}
	N-=bnd;
    pad = bnd/2;
	
	timeSpeed = 1.0f;
	diff = 0.0f;
	visc = 0.0f;
	force = 1000.0f;
	source = 100.0f;

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

	if (!allocate_data_simd()) return 1;
	clear_data();

	win_x = resolution;
	win_y = resolution;

	//DPRINT("Max index of ZIX: " << ZIX(N, N) << "\nMax size of array:" << (N+2)*(N+2) << ".\n");
	//SIMD::m128_test();
	//SIMD::m128_test2(test);
	//return 0;
	
	if (profiling)
	{
		dt = timeSpeed / 60.0f;
	}
	else{
		dt = timeSpeed / (float)frameRate;
	}
	if(profiling){
		while(cur_iter < iterations){
			cur_iter++;
			Simulate(optim_mode);
		}
	}
	else{
		init_sfml(resolution);
		GameLoop();
	}
	return 0;
}