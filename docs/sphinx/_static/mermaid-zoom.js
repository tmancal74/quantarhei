document.addEventListener("DOMContentLoaded", function () {
  var active = null;

  function initZoom() {
    document.querySelectorAll("pre.mermaid svg, .mermaid-container svg").forEach(function (svg) {
      if (svg.dataset.zoomInit) return;
      svg.dataset.zoomInit = "true";

      var container = svg.parentElement;
      container.style.position = "relative";
      container.style.overflow = "hidden";
      container.style.cursor = "grab";

      var state = { scale: 1, panX: 0, panY: 0, startX: 0, startY: 0 };
      svg._zoomState = state;

      function applyTransform() {
        svg.style.transform = "translate(" + state.panX + "px," + state.panY + "px) scale(" + state.scale + ")";
        svg.style.transformOrigin = "0 0";
      }

      container.addEventListener("wheel", function (e) {
        e.preventDefault();
        var rect = container.getBoundingClientRect();
        var mx = e.clientX - rect.left;
        var my = e.clientY - rect.top;
        var factor = e.deltaY > 0 ? 0.9 : 1.1;
        var ns = Math.min(Math.max(state.scale * factor, 0.3), 5);
        var r = ns / state.scale;
        state.panX = mx - r * (mx - state.panX);
        state.panY = my - r * (my - state.panY);
        state.scale = ns;
        applyTransform();
      }, { passive: false });

      container.addEventListener("mousedown", function (e) {
        if (e.button !== 0) return;
        active = { svg: svg, state: state, container: container };
        state.startX = e.clientX - state.panX;
        state.startY = e.clientY - state.panY;
        container.style.cursor = "grabbing";
      });

      container.addEventListener("dblclick", function () {
        state.scale = 1;
        state.panX = 0;
        state.panY = 0;
        applyTransform();
      });
    });
  }

  document.addEventListener("mousemove", function (e) {
    if (!active) return;
    var s = active.state;
    s.panX = e.clientX - s.startX;
    s.panY = e.clientY - s.startY;
    active.svg.style.transform = "translate(" + s.panX + "px," + s.panY + "px) scale(" + s.scale + ")";
    active.svg.style.transformOrigin = "0 0";
  });

  document.addEventListener("mouseup", function () {
    if (!active) return;
    active.container.style.cursor = "grab";
    active = null;
  });

  var observer = new MutationObserver(function (mutations) {
    for (var i = 0; i < mutations.length; i++) {
      if (mutations[i].addedNodes.length) { initZoom(); return; }
    }
  });
  observer.observe(document.querySelector("main") || document.body, { childList: true, subtree: true });
  setTimeout(initZoom, 500);
});
