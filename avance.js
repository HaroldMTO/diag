function avance(val)
{
	var i,j,k,step0,steps,imgs,titles,maps,s;

	maps = document.getElementsByName("map");
	imgs = document.getElementsByName("fig");
	titles = document.getElementsByName("title");
	steps = document.getElementsByName("step");

	if (val == 0) {
		steps[0].options[0].selected = true;
		steps[0].options[0].selectedIndex = 0;
		return;
	}

	step0 = steps[0];
	i = step0.selectedIndex;

	s = step0.options[i].text;

	switch (val) {
	case -2:
		j = 0;
		break;
	case 2:
		j = step0.length-1;
		break;
	default:
		if (i+val < 0) j = 0;
		else if (i+val >= step0.length) j = step0.length-1;
		else j = i+val;
		break;
	}

	if (j == i) return;

	for (k=0;k<imgs.length;k++) imgs[k].src = maps[k].options[j].text;

	for (k=0;k<titles.length;k++) titles[k].innerHTML = step0.options[j].text;

	steps[0].options[i].selected = false;
	steps[0].options[j].selected = true;
	steps[0].selectedIndex = j;
}

