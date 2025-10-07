document.addEventListener("DOMContentLoaded", function () {
  // Get all navigation toggles
  const navToggles = document.querySelectorAll(".nav-toggle");

  navToggles.forEach(function (toggle) {
    toggle.addEventListener("click", function () {
      // Get the associated nav-links
      const navLinks = this.nextElementSibling;
      const arrow = this.querySelector(".arrow");

      // Toggle expanded class (not active)
      this.classList.toggle("expanded");
      navLinks.classList.toggle("expanded");

      // Rotate arrow - CSS should handle this, but keeping for backwards compatibility
      if (this.classList.contains("expanded")) {
        arrow.style.transform = "rotate(90deg)";
      } else {
        arrow.style.transform = "rotate(0deg)";
      }
    });
  });

  // Auto-expand section containing current page
  const currentLinks = document.querySelectorAll(".nav-links a.current");

  currentLinks.forEach(function (currentLink) {
    // Find parent section and expand it
    const parentSection = currentLink.closest(".nav-section");
    if (parentSection) {
      const toggle = parentSection.querySelector(".nav-toggle");
      const navLinksContainer = parentSection.querySelector(".nav-links");
      const arrow = toggle.querySelector(".arrow");

      // Expand the section using 'expanded' class to match CSS
      toggle.classList.add("expanded");
      navLinksContainer.classList.add("expanded");
      arrow.style.transform = "rotate(90deg)";
    }
  });
});
