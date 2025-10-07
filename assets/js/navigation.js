document.addEventListener("DOMContentLoaded", function () {
  // Get all navigation toggles
  const navToggles = document.querySelectorAll(".nav-toggle");

  navToggles.forEach(function (toggle) {
    toggle.addEventListener("click", function () {
      // Get the associated nav-links
      const navLinks = this.nextElementSibling;
      const arrow = this.querySelector(".arrow");

      // Toggle active class
      this.classList.toggle("active");
      navLinks.classList.toggle("active");

      // Rotate arrow
      if (this.classList.contains("active")) {
        arrow.style.transform = "rotate(90deg)";
      } else {
        arrow.style.transform = "rotate(0deg)";
      }
    });
  });

  // Auto-expand section containing current page
  const currentPath = window.location.pathname;
  const navLinks = document.querySelectorAll(".nav-links a");

  navLinks.forEach(function (link) {
    if (
      link.getAttribute("href") === currentPath ||
      currentPath.includes(link.getAttribute("href"))
    ) {
      // Add active class to current link
      link.classList.add("current");

      // Find parent section and expand it
      const parentSection = link.closest(".nav-section");
      if (parentSection) {
        const toggle = parentSection.querySelector(".nav-toggle");
        const navLinksContainer = parentSection.querySelector(".nav-links");
        const arrow = toggle.querySelector(".arrow");

        // Expand the section
        toggle.classList.add("active");
        navLinksContainer.classList.add("active");
        arrow.style.transform = "rotate(90deg)";
      }
    }
  });
});
